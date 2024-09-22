# Struct for holding domain information
class Domain:
    def __init__(self, uid, start, end, seq, stm):
        self.uid = uid
        self.start = start
        self.end = end
        self.seq = seq
        self.stm = stm

    def __str__(self):
        return f"{self.uid}: {self.start}-{self.end}"

    def __repr__(self):
        return str(self)