using Pkg
using PkgTemplates
#Hi
# Defining package with desired plugins
t = Template(user = "RipZyzz2011", plugins = [GitHubActions(), Codecov()])

generate("my_example_pkg",t)