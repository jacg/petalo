module Petalo

double(x) = 2x

using ResumableFunctions

@resumable function onetwothree()
    @yield 1
    @yield 2
    @yield 3
end

end
