class IdentityHolder:
    def __init__(self, value):
        self.value = value
    def __hash__(self):
        return id(self.value)
    def __eq__(self, other):
        return self.value is other.value
    def order(self, other):
        return (self, other) if id(self.value) < id(other.value) else (other, self)
        