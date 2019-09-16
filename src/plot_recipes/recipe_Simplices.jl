using StateSpaceReconstruction

@recipe function f(s::T) where {T <: AbstractSimplex}
    seriestype := :path
    label --> ""
    lc --> :black
    splitaxes(connectvertices(s))
end

@recipe function f(simplices::Vararg{T,N}) where {T <: AbstractSimplex, N}
    seriestype := :path
    legend --> false

    for simplex in simplices
        @series begin
            label --> ""
            lc --> :black
            lw --> 0.8
            splitaxes(connectvertices(simplex))
        end
    end

end


@recipe function f(simplices::AbstractVector{T}) where {T <: AbstractSimplex}
    seriestype := :path
    legend --> false

    for simplex in simplices
        @series begin
            label --> ""
            lc --> :black
            lw --> 0.8
            splitaxes(connectvertices(simplex))
        end
    end

end

export f
