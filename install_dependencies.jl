if VERSION >= v"0.7.0-"
    # Adding Pkg in test/REQUIRE would be an error in 0.6.  Using
    # Project.toml still has some gotchas.  So:
    Pkg = Base.require(Base.PkgId(Base.UUID(0x44cfe95a1eb252eab672e2afdf69b78f), "Pkg"))
end

# Let PyCall.jl use Python interpreter from Conda.jl
# See: https://github.com/JuliaPy/PyCall.jl
ENV["PYTHON"] = ""
Pkg.build("PyCall")

using Conda
Conda.add("scipy")
