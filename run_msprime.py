import click
import msprime
import numpy as np
import tskit

# FIXME: we need to write our own ms headers!
# NOTE: there must be a way to have less repetition for each command?


def infinitely_many_sites_mutation(ts, theta, length, Ne, seed):
    ts = msprime.sim_mutations(
        ts,
        rate=theta / 4.0 / Ne / length,
        model=msprime.BinaryMutationModel(),
        discrete_genome=False,  # needed for "infinite-sites"
        random_seed=seed,
    )
    return ts


def sweep_model(
    nsam, rho, length, Ne, alpha, sweep_age, start_frequency, position, seed
):
    sweep = msprime.SweepGenicSelection(
        position=position,
        start_frequency=start_frequency,
        end_frequency=1 - 1 / 2 / Ne,
        # Because msprime uses s that is 1/2 that of discoal,
        # we only divide by Ne here
        s=alpha / Ne,
        dt=1 / 100 / Ne,
    )
    ts = msprime.sim_ancestry(
        samples=nsam // 2,  # Ack...
        population_size=Ne,
        model=[
            msprime.StandardCoalescent(duration=sweep_age * Ne),
            sweep,
            msprime.StandardCoalescent(),
        ],
        sequence_length=length,
        recombination_rate=rho / 4 / Ne / length,
        random_seed=seed,
        discrete_genome=True,
    )
    return ts


@click.group()
def cli():
    pass


@click.command(name="neutral")
@click.option("--nsam", default=None, type=int)
@click.option("--nreps", default=None, type=int)
@click.option("--theta", default=None, type=float)
@click.option("--rho", default=None, type=float)
@click.option("--length", default=None, type=float)
@click.option("--prefix", default=None, type=str)
@click.option("--seed", default=None, type=int, help="seed")
def neutral(nsam, nreps, theta, rho, length, prefix, seed):
    write_header = True
    open_file = "w"
    for rep, ts in enumerate(
        msprime.sim_ancestry(
            samples=nsam // 2,  # Ack...
            num_replicates=nreps,
            recombination_rate=rho / 4.0 / length,
            sequence_length=length,
            discrete_genome=True,
            random_seed=seed,
        )
    ):
        ts = infinitely_many_sites_mutation(ts, theta, length, 1.0, seed)
        with open(f"{prefix}.ms", open_file) as ofile:
            tskit.write_ms(
                ts, ofile, precision=15, write_header=write_header, num_replicates=nreps
            )
        write_header = False
        open_file = "a"


@click.command(name="soft")
@click.option("--nsam", default=None, type=int, help="Sample size.")
@click.option("--nreps", default=None, type=int, help="Number of replicates.")
@click.option("--theta", default=None, type=float, help="theta = 4*Ne*mu")
@click.option("--rho", default=None, type=float, help="rho = 4*Ne*r")
@click.option("--length", default=None, type=float, help="Genome length")
@click.option("--prefix", default=None, type=str, help="Output file name.")
@click.option("--Ne", default=None, type=float, help="Effective population size")
@click.option("--alpha", default=None, type=float, help="2Ne*s")
@click.option(
    "--freq",
    default=None,
    type=float,
    help="Starting frequency of sweep.",
)
@click.option("--position", default=None, type=float, help="Position of sweep.")
@click.option("--seed", default=None, type=int, help="seed")
def soft_sweep(
    nsam,
    nreps,
    theta,
    rho,
    length,
    prefix,
    ne,
    alpha,
    freq,
    position,
    seed,
):
    Ne = ne

    write_header = True
    open_file = "w"

    used_seeds = {}
    generator = np.random.Generator(np.random.MT19937(seed))
    for _ in range(nreps):
        sweep_age = 0.0

        seed = generator.integers(0, np.iinfo(np.uint32).max)

        while seed in used_seeds:
            seed = generator.randint(0, hi=np.iinfo(np.uint32).max)

        used_seeds[seed] = 1

        ts = sweep_model(
            nsam,
            rho,
            length,
            Ne,
            alpha,
            sweep_age,
            freq,
            position,
            seed,
        )

        ts = infinitely_many_sites_mutation(ts, theta, length, Ne, seed)
        with open(f"{prefix}.ms", open_file) as ofile:
            tskit.write_ms(
                ts, ofile, precision=15, write_header=write_header, num_replicates=nreps
            )
        write_header = False
        open_file = "a"


cli.add_command(soft_sweep)
cli.add_command(neutral)

if __name__ == "__main__":
    cli()
