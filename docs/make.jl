using Documenter, VLQuantitativeFinancePackage

push!(LOAD_PATH,"../src/")

makedocs(
    sitename="VLQuantitativeFinancePackage.jl",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [VLQuantitativeFinancePackage],
    pages = [
        "Home" => "index.md",
        "Treasury securities" => "fixed.md",
        "Equity securities" => "equity.md",
        "Derivative securities" => "derivatives.md",
        "Portfolio management" => "portfolio.md",
    ]
)

deploydocs(
    repo = "github.com/varnerlab/VLQuantitativeFinancePackage.jl.git", branch = "gh-pages", target = "build"
)