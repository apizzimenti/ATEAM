
class Vertex:
    def __init__(self, at, spin, index):
        self.at = at
        self.spin = spin
        self.index = index

    def __eq__(self, other): return self.at == other.at

    def __hash__(self): return hash(self.at)

    def __str__(self):
        return str(self.at) + f"\t\t spin = {self.spin}\t\t index = {self.index}"


class Edge:
    def __init__(self, at, spin, index):
        self.at = at
        self.spin = spin
        self.index = index

    def asIndices(self): return self.at[0].index, self.at[1].index

    def __eq__(self, other): return self.at == other.at

    def __hash__(self): return hash(self.at)

    def __str__(self):
        return f"({str(self.at[0].at)},{self.at[1].at})\t\t spin = {self.spin}\t\t index = {self.index}"
