from sfbox_api import Frame, Lat, Mol, Mon, Sys

if __name__ == "__main__":
    alpha = 0.1
    phi_bulk = 1e-3

    lat = Lat(**{"geometry": "flat", "n_layers": 100})

    mols = []
    mols.append(Mol(**{"name": "water", "composition": "(W)1", "freedom": "solvent"}))
    mols.append(
        Mol(
            **{
                "name": "pol",
                "composition": "(X)1(A)98(G)1",
                "freedom": "restricted",
                "theta": 5.0,
            }
        )
    )
    mols.append(
        Mol(
            **{
                "name": "ionp",
                "composition": "(P)1",
                "freedom": "free",
                "phibulk": phi_bulk,
            }
        )
    )
    mols.append(
        Mol(**{"name": "ionm", "composition": "(M)1", "freedom": "neutralizer"})
    )

    mons = []
    mons.append(Mon(**{"name": "W", "freedom": "free"}))
    mons.append(
        Mon(**{"name": "X", "freedom": "pinned", "pinned_range": "1", "valence": alpha})
    )
    mons.append(Mon(**{"name": "A", "freedom": "free", "valence": alpha}))
    mons.append(Mon(**{"name": "G", "freedom": "free", "valence": alpha}))

    mons.append(Mon(**{"name": "P", "freedom": "free", "valence": 1.0}))
    mons.append(Mon(**{"name": "M", "freedom": "free", "valence": -1.0}))

    sys = Sys()

    frame = Frame(lat, sys, mols, mons)

    frame.run()

    print(frame.profile_labels)
    # print(frame.profile["water"])
    # print(frame.profile["name"])