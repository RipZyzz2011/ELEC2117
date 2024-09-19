using Pkg
using PkgTemplates

t = Template(user = "RipZyzz2011", plugins = [GitHubActions(), Codecov()])

generate("LinearODESolver", t)