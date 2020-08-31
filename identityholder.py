class IdentityHolder:
    def __init__(self, value, tag=None):
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
        return (self, other) if id(self.value) < id(other.value) else (other, self)
    def isbefore(self, other):
        return id(self.value) < id(other.value)
        