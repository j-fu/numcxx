import matplotlib


def plotGeometry(plt,geom):
    points=geom.get_points()
    bfaces=geom.get_bfaces()
    bfaceregions=geom.get_bfaceregions()
    regionpoints=geom.get_regionpoints()
    regionnumbers=geom.get_regionnumbers()

    colortable=[ 
        (0,0,0),
        (1,0,0),
        (0,1,0),
        (0,0,1),
        (1,1,0),
        (0,1,1),
        (1,0,1)]

    for ipoint in range(regionpoints.shape[0]):
        plt.plot([regionpoints[ipoint][0]],[regionpoints[ipoint][1]],
                 'o-',color=colortable[regionnumbers[ipoint]%len(colortable)])

    for ipoint in range(points.shape[0]):
        plt.plot([points[ipoint][0]],[points[ipoint][1]],
                 'o-',color=colortable[0])


    for ibface in range(bfaces.shape[0]):
        p=[bfaces[ibface][0],bfaces[ibface][1]]
        plt.plot([points[p[0]][0],points[p[1]][0]],
                 [points[p[0]][1],points[p[1]][1]],
                 '-',color=colortable[bfaceregions[ibface]%len(colortable)])

def plotGrid(plt,grid):
    points=grid.get_points()
    cells=grid.get_cells()
    bfaces=grid.get_bfaces()
    bfaceregions=grid.get_bfaceregions()

    colortable=[ 
        (0,0,0),
        (1,0,0),
        (0,1,0),
        (0,0,1),
        (1,1,0),
        (0,1,1),
        (1,0,1)]

    for icell in range(grid.ncells()):
        p=[cells[icell][0],cells[icell][1],cells[icell][2]]
        plt.plot([points[p[0]][0],points[p[1]][0],points[p[2]][0],points[p[0]][0]],
                 [points[p[0]][1],points[p[1]][1],points[p[2]][1],points[p[0]][1]],
                 'go-')

    for ibface in range(grid.nbfaces()):
        p=[bfaces[ibface][0],bfaces[ibface][1]]
        plt.plot([points[p[0]][0],points[p[1]][0]],
                 [points[p[0]][1],points[p[1]][1]],
                 'o-',color=colortable[bfaceregions[ibface]%len(colortable)])




def triangulation(grid):
    xy=grid.get_points()
    x=xy[:,[0]].squeeze()
    y=xy[:,[1]].squeeze()
    triangles=grid.get_cells()
    return matplotlib.tri.Triangulation(x,y,triangles)
