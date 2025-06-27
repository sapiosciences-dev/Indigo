from __future__ import annotations

from typing import Callable, Dict, List, Optional
from uuid import uuid4

from indigo import Indigo, IndigoException, IndigoObject  # type: ignore

MOL_TYPES = ["#02: <molecule>", "#03: <query reaction>", "#12: <RDFMolecule>"]
REAC_TYPES = ["#04: <reaction>", "#05: <query reaction>"]

def sanitize_indigo_query_molecule(query_mol: IndigoObject):
    """
    Replaces pseudo-atoms, R-groups, and wildcard atoms in an Indigo query molecule
    with proper wildcard atoms (*).
    """
    for atom in query_mol.iterateAtoms():
        atoms_to_remove = []
        to_remove = False
        if atom.isPseudoatom() or atom.isRSite():
            to_remove = True
        else:
            atom.removeConstraints("substituents")
            atom.removeConstraints("substituents-as-drawn")
            symbol = atom.symbol()
            if symbol in {'A', 'Q', 'X', 'M', 'AH', 'QH', 'XH', 'MH', 'NOT', 'R', '*'}:
                to_remove = True
            if "[" in symbol and "]" in symbol:
                to_remove = True
        if to_remove:
            atoms_to_remove.append(atom)

        query_mol.removeAtoms([x.index() for x in atoms_to_remove])
    return query_mol

def sanitize_indigo_query_reaction(query_rxn: IndigoObject):
    """
    Remove all query features from reaction. This is to store a simplified version for query-query search.
    """
    indigo: Indigo = query_rxn.session
    cloned = indigo.createQueryReaction()
    for mol in query_rxn.iterateReactants():
        cloned_mol = mol.clone()
        sanitize_indigo_query_molecule(cloned_mol)
        if cloned_mol.countAtoms() > 0:
            cloned.addReactant(cloned_mol)
    for mol in query_rxn.iterateProducts():
        cloned_mol = mol.clone()
        sanitize_indigo_query_molecule(cloned_mol)
        if cloned_mol.countAtoms() > 0:
            cloned.addProduct(cloned_mol)
    if cloned.countProducts() == 0 and cloned.countReactants() == 0:
        raise ValueError("At least one reactant or product must contains a non-query atom in the reaction template.")
    return indigo.loadReaction(cloned.rxnfile())

# pylint: disable=unused-argument
def skip_errors(instance: IndigoRecord, err: BaseException) -> None:
    """
    Empty handler to skip errors
    """


def check_error(instance: IndigoRecord, error: BaseException) -> None:
    if instance.error_handler:
        instance.error_handler(instance, error)
    else:
        raise error


class WithElasticResponse:
    def __set__(self, instance: IndigoRecord, value: Dict):
        el_src = value["_source"]
        for arg, val in el_src.items():
            setattr(instance, arg, val)
        setattr(instance, "_sort", value.get("sort"))


class WithIndigoObject:
    def __set__(  # pylint: disable=too-many-branches
        self, instance: IndigoRecord, value: IndigoObject
    ) -> None:
        try:
            value_dup = value.clone()
        except IndigoException as err_:
            check_error(instance, err_)
            # Can't create IndigoObject
            return

        value_dup.aromatize()
        fingerprints = (
            "sim",
            "sub",
        )

        for f_print in fingerprints:
            try:
                setattr(instance, f"{f_print}_fingerprint", [])
                setattr(instance, f"{f_print}_fingerprint_len", 0)

                fp_list = value_dup.fingerprint(f_print).oneBitsList()
                if fp_list:
                    fp_ = [int(feature) for feature in fp_list.split(" ")]
                    setattr(instance, f"{f_print}_fingerprint", fp_)
                    setattr(instance, f"{f_print}_fingerprint_len", len(fp_))
            except ValueError as err_:
                check_error(instance, err_)
            except IndigoException as err_:
                check_error(instance, err_)

        try:
            cmf = " ".join(map(str, list(value_dup.serialize())))
            setattr(instance, "cmf", cmf)
        except IndigoException as err_:
            setattr(instance, "cmf", "")
            check_error(instance, err_)

        try:
            setattr(instance, "name", value_dup.name())
        except IndigoException as err_:
            setattr(instance, "name", "")
            check_error(instance, err_)

        try:
            internal_type = value_dup.dbgInternalType()
            if internal_type in MOL_TYPES:
                hash_ = [
                    component.clone().hash()
                    for component in value_dup.iterateComponents()
                ]
                setattr(instance, "hash", sorted(set(hash_)))
            elif internal_type in REAC_TYPES:
                setattr(instance, "hash", [value_dup.hash()])
        except IndigoException as err_:
            check_error(instance, err_)

        try:
            if value_dup.checkBadValence():
                setattr(instance, "has_error", 1)
            else:
                setattr(instance, "has_error", 0)
        except IndigoException as err_:
            check_error(instance, err_)


class IndigoRecord:
    """
    Base class for IndigoObject representation.
    This class could not be instantiated directly, use one of the following
    subclasses:
        - IndigoRecordMolecule
        - IndigoRecordReaction
        - IndigoRecordReactionTemplate
    """

    cmf: Optional[str] = None
    name: Optional[str] = None
    rawData: Optional[str] = None
    sim_fingerprint: Optional[List[str]] = None
    sub_fingerprint: Optional[List[str]] = None
    indigo_object = WithIndigoObject()
    elastic_response = WithElasticResponse()
    record_id: Optional[str] = None
    error_handler: Optional[Callable[[object, BaseException], None]] = None
    # Sort for page cursor
    _sort: list

    @property
    def sort(self) -> list:
        return self._sort

    def __init__(self, **kwargs) -> None:
        """
        Constructor accepts only keyword arguments
        :param indigo_object: — create indigo record from IndigoObject
        :type indigo_object: IndigoObject
        :param name: — add name. Rewrites name from IndigoObject
        :type name: str
        :param sim_fingerprint: similarity fingerprint (sim)
        :type sim_fingerprint: List[int]
        :param sub_fingerprint: similarity fingerprint (sub)
        :type sub_fingerprint: List[int]
        :param error_handler: lambda for catching exceptions
        :type error_handler: Optional[Callable[[object, BaseException], None]]
        :param skip_errors: if True, all errors will be skipped,
                            no error_handler is required
        :type skip_errors: bool
        """

        # First check if skip_errors flag passed
        # If no flag passed add error_handler function from arguments
        if kwargs.get("skip_errors", False):
            self.error_handler = skip_errors  # type: ignore
        else:
            self.error_handler = kwargs.get("error_handler", None)

        self.record_id = uuid4().hex
        for arg, val in kwargs.items():
            setattr(self, arg, val)

    def as_dict(self) -> Dict:
        # Add system fields here to exclude from indexing
        filtered_fields = {"error_handler", "skip_errors"}
        return {
            key: value
            for key, value in self.__dict__.items()
            if key not in filtered_fields
        }

    def as_indigo_object(self, session: Indigo):
        if self.cmf:
            return session.deserialize(list(map(int, self.cmf.split(" "))))  # type: ignore
        raise ValueError("Unexpected cmf value")


class IndigoRecordMolecule(IndigoRecord):
    pass


class IndigoRecordReaction(IndigoRecord):
    pass

class IndigoRecordReactionTemplate(IndigoRecord):

    def __init__(self, **kwargs) -> None:
        """
        Convert the query reaction into a cleaned up version with wildcards other query features removed.
        Then store the original in the cartridge as regular text RXN.
        We can reverse the match if the user later provides a concrete reaction for better accuracy.
        """
        if "indigo_object" in kwargs:
            original_query_reaction: IndigoObject = kwargs["indigo_object"]
            cleaned_reaction: IndigoObject = sanitize_indigo_query_reaction(original_query_reaction)
            kwargs["indigo_object"] = cleaned_reaction
            kwargs["original_qrxn"] = original_query_reaction.rxnfile()
        super().__init__(**kwargs)



def as_iob(indigo_record: IndigoRecord, session: Indigo) -> IndigoObject:
    """Function extracts IndigoObject from IndigoRecord
    Short alias to IndigoRecord.as_indigo_object
    :param indigo_record:
    :param session:
    :return:
    """
    return indigo_record.as_indigo_object(session)
