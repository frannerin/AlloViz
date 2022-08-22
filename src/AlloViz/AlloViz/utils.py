"""
Module containing helper functions and variables used by many other modules

"""

from .info import wrappers

pkgsl = list(wrappers.keys())
"""List of the correct names of available network construction methods

Imported from :data:`AlloViz.AlloViz.info.wrappers`"""


def pkgname(pkg, fail=True):
    r"""Return the case-sensitive, correct name of an AlloViz network construction method

    Check if the passed string matches any of the available AlloViz network construction
    methods detailed in :data:`AlloViz.AlloViz.info.wrappers` and return the correctly
    written AlloViz accession name, else raise an Exception or return False.

    Parameters
    ----------
    pkg : str
        Name for which to retrieve the correct AlloViz accession name.
    fail : bool, default=True
        Raise an Exception if an AlloViz accession name cannot be retrieved or simply
        return False.
    """
    lowpkgsl = [p.lower() for p in pkgsl]
    pkgname = pkgsl[lowpkgsl.index(pkg.lower())] if pkg.lower() in lowpkgsl else False

    if not pkgname and fail:
        raise Exception(
            f"{pkg} isn't a valid name of an AlloViz network construction method."
        )

    return pkgname


#: List of the available network metrics that can be calculated in analysis
#:
#: To be used when the metrics="all" parameter is used in a function
metricsl = ["cfb", "btw"]

#: List of the available network filterings that can be performed for analysis
#:
#: To be used when the filterby="all" parameter is used in a function
filteringsl = ["All", "GetContacts_edges", "No_Sequence_Neighbors", "GPCR_Interhelix"]


def make_list(obj, if_all, apply=lambda x: x):
    r"""Process the passed object and return a list of strings
    
    Check if the object is a string, a list, or "all" and return the a list with the
    appropriate values (if it is the case, after applying the `apply` function to them).
    Make the string a list of length 1, return the unedited list, or return the passed 
    `if_all` object if `obj` is '"all"'. If it is the case, the individual string(s) of
    the returned list are processed with the passed `apply` function (default does
    nothing).

    Parameters
    ----------
    obj : str or list
        String or list to process.
    if_all
        Object to return if `obj == "all"`.
    apply : func, optional
        Function to apply to the list's strings before returning. Default does nothing.
    """
    return if_all if obj == "all" else [apply(o) for o in obj] if isinstance(obj, list) else [apply(obj)]


def rgetattr(obj, *attrs):
    r"""Recursive version of the built-in getattr

    It recursively checks if the successive strings passed are attributes of the object
    or the object's attribute, or the object's attribute's attribute... to finally return
    the final attribute or else return False.

    Parameters
    ----------
    obj
        Object in which to check if the first attribute passed exists and retrieve it.
    attrs : str
        Strings to recursively use to retrieve attributes.

    See also
    --------
    AlloViz.AlloViz.utils.rhasattr

    Examples
    --------
    >>> import AlloViz
    >>> rgetattr = AlloViz.AlloViz.utils.rgetattr
    >>> rgetattr
    <function AlloViz.AlloViz.utils.rgetattr(obj, *attrs)>
    >>> rgetattr(AlloViz, "AlloViz", "utils", "rgetattr")
    <function AlloViz.AlloViz.utils.rgetattr(obj, *attrs)>
    >>> rgetattr(AlloViz.AlloViz.utils, "rgetattr")
    <function AlloViz.AlloViz.utils.rgetattr(obj, *attrs)>
    >>> rgetattr(AlloViz, "AlloViz", "utils", "is_attr")
    False
    """
    # Make a list of all the strings passed
    attl = list(attrs)
    # While there are still strings in the list, try to retrieve the attribute from 'obj' and rewrite 'obj' with it to keep recursing, else return False if it doesn't exist
    while len(attl) >= 1:
        attr = attl.pop(0)
        try:
            obj = getattr(obj, attr)
        except AttributeError:
            return False
    return obj


def rhasattr(obj, *attrs):
    r"""Recursive version of the built-in hasattr

    It recursively checks if the successive strings passed are attributes of the object
    or the object's attribute, or the object's attribute's attribute... It exploits
    :func:`~AlloViz.AlloViz.utils.rgetattr` and the fact that it already returns False if
    any of the strings passed for the recursive search doesn't exist as attribute.

    Parameters
    ----------
    obj
        Object in which to check if the first attribute passed exists.
    attrs : str
        Strings to recursively use to retrieve attributes.

    Examples
    --------
    >>> import AlloViz
    >>> rgetattr = AlloViz.AlloViz.utils.rgetattr
    >>> rgetattr(AlloViz, "AlloViz", "utils", "rgetattr")
    True
    >>> rgetattr(AlloViz.AlloViz.utils, "rgetattr")
    True
    >>> rgetattr(AlloViz, "AlloViz", "utils", "is_attr")
    False
    """
    # Return True if the output of rgetattr is anything but a boolean (might lead to error if the checked attribute is a boolean itself)
    # Else, if it is a boolean it is assumed to be a False return caused by an unexisting attribute, and thus False is returned
    result = rgetattr(obj, *attrs)
    return True if not isinstance(result, bool) else False


class dummypool:
    r"""Class to mimic a process Pool when only using 1 core (synchronous)

    This class aims to be able to be used with the same syntax as a multiprocess(ing)
    Pool managed through ``apply_async``. Instead of running the tasks sent with the
    method asynchronously, if the pool is an instance of this class they will run
    synchronously in each call to the method.

    Notes
    -----
    Using a multiprocess(ing) Pool initialized with 1 core would have the same effect in
    terms of resource consumption, but this way the tasks are run in the main namespace
    instead of on a pickled copy, which is useful for debugging, i.e., getting the
    stdout and stderr immediately.
    """

    def apply_async(self, function, args=[], callback=None):
        r"""Execute the function with the passed args

        It executes the function with the passed args synchronously, and optionally
        the specified callback function to the resulting output as well.

        Parameters
        ----------
        function : func
        args : list
        callback : func, optional
        """
        if callback is not None:
            return callback(function(*args))
        else:
            return function(*args)

    def close(self):
        """Empty function"""
        pass

    def join(self):
        """Empty function"""
        pass


#: Pool variable to share among modules
#:
#: Defining a `pool` variable inside this module allows for other modules to modify it
#: and share it between modules, even when pickling due to the use of a multiprocess.Pool
pool = dummypool()


def get_pool():
    r"""Function to retrieve shared pool variable

    This function retrieves this module's pool variable from the main namespace and
    returns it to use it in whatever namespace it is called from.
    """
    global pool
    return pool


# def capitalize(string):
#     return string[0].upper() + string[1:]

# def norm(normalize):
#     return "norm" if normalize else "no_norm"


# def get_intercontacts(indexl):
#     get_resnum = lambda res: int(res.rsplit(":")[-1])
#     return [idx for idx in indexl if abs(get_resnum(idx[0]) - get_resnum(idx[1]) ) >= 4]


# pdict = {}

# def update_pdict(name, p):
#     global pdict
#     pdict[name] = p
#     print("adding to pdict", pdict)
#     p.start()

# def get_p(name):
#     global pdict
#     print("getting from pdict", pdict)
#     p = pdict[name]
#     print(p)
#     p.get()
#     print(p)
