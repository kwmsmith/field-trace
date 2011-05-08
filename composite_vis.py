import pylab as pl
def visualize(
        arr,
        mesh=None,
        crit_pts=None,
        surf_network=None,
        len_lim=None,
        cmap=None,
        ncontours=None,
        new_fig=True,
        save_fig=None,
        cpt_size=20,
        fig_title=None,
        exts=('.eps', '.png', '.pdf')
        ):

    if new_fig:
        fig = pl.figure()
    # pl.imshow(arr, interpolation='nearest', cmap='jet')
    pl.imshow(arr, cmap=cmap, interpolation='nearest')
    pl.colorbar(pad=0.0, shrink=0.9)
    if ncontours is not None:
        pl.contour(arr, ncontours, linewidths=2, cmap=pl.cm.hot)
    nx, ny = arr.shape
    if mesh is not None:
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                other_pts = mesh._g[i,j]
                for x,y in other_pts:
                    pl.plot([j,y], [i,x], 'k--')
    if surf_network is not None:
        for node in surf_network:
            node_x, node_y = node
            nbrs = surf_network[node]
            for nbr in nbrs:
                nbr_x, nbr_y = nbr
                if len_lim and abs(node_y-nbr_y) < nx /2 and abs(node_x-nbr_x) < ny /2:
                    pl.plot([node_y, nbr_y], [node_x, nbr_x], 'k--', linewidth=2)
    if crit_pts is not None:
        pits, passes, peaks = crit_pts.pits, crit_pts.passes, crit_pts.peaks
        X = [_[0] for _ in pits]
        Y = [_[1] for _ in pits]
        pl.scatter(Y, X, marker='o', c='b', s=cpt_size)
        X = [_[0] for _ in peaks]
        Y = [_[1] for _ in peaks]
        pl.scatter(Y, X, marker='s', c='r', s=cpt_size)
        X = [_[0] for _ in passes]
        Y = [_[1] for _ in passes]
        pl.scatter(Y, X, marker='d', c='k', s=cpt_size)
    pl.xlim(0, nx)
    pl.ylim(0, ny)
    if fig_title:
        pl.title(fig_title)
    if new_fig:
        pl.figure(fig.number)
    else:
        pl.figure(pl.gcf().number)
    if save_fig:
        for ext in exts:
            pl.savefig(save_fig+ext)
