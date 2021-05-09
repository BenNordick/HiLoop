class IdentityHolder:
    """
    Wraps an object so that it can be quickly compared by object identity rather than value.
    Also provides a total order based on object address.

    This is frequently useful for sets representing cycles: set value equality is O(n) but object identity equality is O(1).
    Different cycles can have the same nodes in a different order; it is important to distinguish these objects, but useful
    to represent them as sets for purposes of fast intersection.
    """
    def __init__(self, value, tag=None):
        """Wrap a value, optionally tagged with extra information."""
        self.value = value
        self.tag = tag
    def __hash__(self):
        return id(self.value)
    def __eq__(self, other):
        return self.value is other.value
    def __str__(self):
        return 'IH`' + str(self.value) + '`'
    def __repr__(self):
        return 'IdentityHolder(' + repr(self.value) + ', ' + repr(self.tag) + ')'
    def order(self, other):
        """Return a tuple of this holder and the other such that the first item is before the second in a total order."""
        return (self, other) if self.isbefore(other) else (other, self)
    def isbefore(self, other):
        """Determine whether this holder is before the other in a total order."""
        return id(self.value) < id(other.value)
        