using Pkg
Pkg.activate("docs")
Pkg.develop(PackageSpec(path="."))
Pkg.instantiate()
