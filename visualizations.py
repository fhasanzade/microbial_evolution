import gzip, pickle, sys
import matplotlib.pyplot as plt


def main():
    fname = sys.argv[1]

    fh = gzip.open(fname, "rb")
    params = pickle.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T, V, H, I = [], [], [], []
    for epoch in range(params["epochs"]):
        viruses = pickle.load(fh)
        hosts = pickle.load(fh)
        interactions = pickle.load(fh)
        T.append(epoch)
        V.append(len(viruses))
        H.append(len(hosts))
        I.append(len(interactions))
    fh.close()

    figure = plt.figure(figsize=(8, 6), dpi=500)
    size = 6
    plt.rcParams["axes.titlesize"] = size
    plt.rcParams["axes.labelsize"] = size
    plt.rcParams["xtick.labelsize"] = size
    plt.rcParams["ytick.labelsize"] = size
    plt.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(3, 1, 1)
    axes.set_ylabel("# of viruses")
    axes.plot(T, V, "r-", alpha = 0.6)

    axes = figure.add_subplot(3, 1, 2)
    axes.set_ylabel("# of hosts")
    axes.plot(T, H, "b-", alpha = 0.6)

    axes = figure.add_subplot(3, 1, 3)
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("# of virus/host interactions")
    axes.plot(T, I, "g-", alpha = 0.6)

    plt.savefig("%s" % (fname.replace("pkl", "pdf")), format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
