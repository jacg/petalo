module Petalo

export synchronize

using ResumableFunctions

mutable struct Synchronizer
    iterable
    key
    data
    count
    iterator_state

    function Synchronizer(iterable, key)
        s = new(iterable, key, nothing, nothing, nothing)
        s.data, s.iterator_state = iterate(iterable)
        s.count = s.key(s.data)
        s
    end
end


next(s::Synchronizer) = next(s, iterate(s.iterable, s.iterator_state))
next( ::Synchronizer, ::Nothing) = false
function next(s::Synchronizer, (data, state))
    s.data, s.iterator_state = data, state
    s.count = s.key(s.data)
    return true
end



getcount(x) = getproperty(x, :count)
getdata(x)  = getproperty(x, :data)
mincount(sources) = minimum(getcount, sources)
maxcount(sources) = maximum(getcount, sources)

@resumable function synchronize(iterables, keys)
    sources = map(args -> Synchronizer(args...), zip(iterables, keys))
    while true
        while mincount(sources) < (highest_count = maxcount(sources))
            for source in sources
                if source.count < highest_count
                    next(source)
                end
            end
        end
        @yield highest_count, map(getdata, sources)
        for source in sources
            if ! next(source) return end
        end
    end
end

end
