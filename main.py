from dataclasses import dataclass


@dataclass
class Test:
    number: int
    string: str


print(Test(1, "1"))
print(list(filter(lambda test: test.number == 2, [Test(1, "1"), Test(2, "2"), Test(3, "3")])))
