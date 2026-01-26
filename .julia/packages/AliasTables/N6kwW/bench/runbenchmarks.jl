using RegressionTests, Chairmarks

load_time = @elapsed using AliasTables
@track load_time

x = rand(10)
ttf_construction = @elapsed at = AliasTable(x)
@track ttf_construction

ttf_sample = @elapsed rand(at)
@track ttf_sample

ttf_sample_multi = @elapsed rand(at, 30)
@track ttf_sample_multi

for n in [1, 10, 100, 1000]
    @track @b rand(n) AliasTable seconds=.02
    @track @b AliasTable(rand(n)) rand seconds=.02
    @track @b AliasTable(rand(n)) rand(_, 100) seconds=.02
end

for n in [3, 30, 300]
    @track @b AliasTable(rand(n)) hash seconds=.01
    @track @b AliasTable(rand(n)), AliasTable(rand(n)) (==)(_...) seconds=.01
    @track @b (x = AliasTable(rand(n)); (x, deepcopy(x))) (==)(_...) seconds=.01
end
