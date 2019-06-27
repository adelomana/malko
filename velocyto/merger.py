import loompy

files=[
    '/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/20765_M397_control/velocyto/20765_M397_control.loom',
    '/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/21729_day3/velocyto/21729_day3.loom',
    '/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/21766_day6/velocyto/21766_day6.loom',
    '/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/22077_day13/velocyto/22077_day13.loom',
    '/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/22152_day17/velocyto/22152_day17.loom',
    '/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/18324_Yapeng_single_cell/velocyto/18324_Yapeng_single_cell.loom'
    ]
output_filename='/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/combined.loom'
loompy.combine(files, output_filename, key="Accession")
