using Pkg
using PkgTemplates

# Defining package with desired plugins
t = Template(user = "RipZyzz2011", plugins = [GitHubActions(), Codecov()])

generate("MyStatsPackage",t)