# Comment1
class Module1(object):
    command_name = "module1"
    targets = [r"module1_target.txt"]
    dependencies = [("module1_dep*.txt", True), ("not_a_dependency.txt", False)]
    configs = ["module1a.conf", "module1b.conf"]
