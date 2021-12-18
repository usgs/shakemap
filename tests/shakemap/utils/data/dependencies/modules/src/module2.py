# Comment1
class Module2(object):
    command_name = "module2"
    targets = [r"products/module2_target.txt"]
    dependencies = [("module1_target.txt", True)]
    configs = ["module2a.conf", "module2b.conf"]
