from ROOT import TFile, TTree
from array import array
from pathlib import Path
import polars as pl
from tqdm import tqdm


workspace_path = Path('/data/a1975/')
solved_path = workspace_path / 'InterpSolver'
particle = '4He'
root_file_path = solved_path / f'solved_{particle}.root'
run_min = 106
run_max = 189


def main():
	# read file
	file_list = []
	for run in range(run_min, run_max+1):
		run_path = solved_path / f'run_{run:04d}_{particle}.parquet'
		if not run_path.exists():
			print(f'Warning: Run {run} not found!')
			continue
		file_list.append(run_path)
	df = pl.read_parquet(file_list)

	# create root file and tree
	root_file = TFile(str(root_file_path), "RECREATE")
	tree = TTree('tree', "Solved data from parquet file")

	# data
	cluster_index = array('l', [0])
	cluster_label = array('l', [0])
	origin_run = array('l', [0])
	origin_event = array('l', [0])
	vertex_x = array('d', [0.0])
	vertex_y = array('d', [0.0])
	vertex_z = array('d', [0.0])
	brho = array('d', [0.0])
	kinetic_energy = array('d', [0.0])
	polar = array('d', [0.0])
	azimuthal = array('d', [0.0])
	reduced_chi2 = array('d', [0.0])
	sigma_vx = array('d', [0.0])
	sigma_vy = array('d', [0.0])
	sigma_vz = array('d', [0.0])
	sigma_brho = array('d', [0.0])
	sigma_kinetic_energy = array('d', [0.0])
	sigma_polar = array('d', [0.0])
	sigma_azimuthal = array('d', [0.0])
	# branches
	tree.Branch('origin_run', origin_run, 'run/I')
	tree.Branch('origin_event', origin_event, 'evt/I')
	tree.Branch('cluster_index', cluster_index, 'cindex/I')
	tree.Branch('cluster_label', cluster_label, 'clabel/I')
	tree.Branch('vertex_x', vertex_x, 'vx/D')
	tree.Branch('vertex_y', vertex_y, 'vy/D')
	tree.Branch('vertex_z', vertex_z, 'vz/D')
	tree.Branch('brho', brho, 'brho/D')
	tree.Branch('kinetic_energy', kinetic_energy, 'ek/D')
	tree.Branch('polar', polar, 'theta/D')
	tree.Branch('azimuthal', azimuthal, 'phi/D')
	tree.Branch('reduced_chi2', reduced_chi2, 'rchi2/D')
	tree.Branch('sigma_vx', sigma_vx, 'svx/D')
	tree.Branch('sigma_vy', sigma_vy, 'svy/D')
	tree.Branch('sigam_vz', sigma_vz, 'svz/D')
	tree.Branch('sigma_brho', sigma_brho, 'sbrho/D')
	tree.Branch('sigma_kinetic_energy', sigma_kinetic_energy, 'sek/D')
	tree.Branch('sigma_polar', sigma_polar, 'stheta/D')
	tree.Branch('sigma_azimuthal', sigma_azimuthal, 'sphi/D')



	for row in tqdm(df.iter_rows(), desc='convert', unit='events'):
		cluster_index[0] = row[1]
		cluster_label[0] = row[2]
		origin_run[0] = row[3]
		origin_event[0] = row[4]
		vertex_x[0] = row[5]
		vertex_y[0] = row[7]
		vertex_z[0] = row[9]
		brho[0] = row[11]
		kinetic_energy[0] = row[13]
		polar[0] = row[15]
		azimuthal[0] = row[17]
		reduced_chi2[0] = row[19]
		sigma_vx[0] = row[6]
		sigma_vy[0] = row[8]
		sigma_vz[0] = row[10]
		sigma_brho[0] = row[12]
		sigma_kinetic_energy[0] = row[14]
		sigma_polar[0] = row[16]
		sigma_azimuthal[0] = row[18]
		tree.Fill()

	# save and close file
	tree.Write()
	root_file.Close()


if __name__ == "__main__":
	main()