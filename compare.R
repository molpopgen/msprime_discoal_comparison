pdf("compare.pdf")
par(mfrow=c(2,2))
mn = scan("msprime_neutral.segsites")
for (i in c(1, 100, 500, 1000))
{
    mspfile = paste("msprime_", i, ".segsites", sep="")
    dscfile = paste("discoal_", i, ".segsites", sep="")
    fp11file = paste("fp11_", i, ".segsites", sep="")

    m = scan(mspfile, quiet=TRUE)
    d = scan(dscfile, quiet=TRUE)
    f = scan(fp11file, quiet=TRUE)

    plot(ecdf(m), main=paste("Initial freq = ", 0.5*i/1e4), xlim=c(2500,5000))
    lines(ecdf(d), col="blue")
    lines(ecdf(f), col="magenta")
    lines(ecdf(mn), col="red")

    print(paste(i, 0.5*i/1e4, mean(m), mean(d), ks.test(m, d)$p.value))
}

