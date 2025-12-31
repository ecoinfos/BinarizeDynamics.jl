using BinarizeDynamics
using Random

try
    seqs = ["AA", "BB"]
    diff_dynamics(seqs, seqs; test=:bootstrap, n_resamples=2)
    println("Success!")
catch e
    open("error_log.txt", "w") do io
        showerror(io, e)
        println(io)
        show(io, stacktrace(catch_backtrace()))
    end
    rethrow(e)
end
