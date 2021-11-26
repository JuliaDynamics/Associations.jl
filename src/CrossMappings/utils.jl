abstract type CrossmapEmbedding end

struct CCMEmbedding end
struct PAIEmbedding end

function crossmapembed(x, d, τ, method::CCMEmbedding)
    τs = 0:-τ:d*-τ+1
    Mₓ = genembed(x, τs); 
    return Mₓ
end

function crossmapembed(x, y, d, τ, method::PAIEmbedding)
    τs = [collect(0:-τ:d*-τ+1); 0]
    js = [repeat([1], d); 2]
    Mₓy = genembed(Dataset(x, y), τs, js); 
    return Mₓy
end