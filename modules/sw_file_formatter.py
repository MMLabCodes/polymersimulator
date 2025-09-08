from modules.sw_basic_functions import *


class DFT_input_generator():
    dft_template = """
    {keywords}
    %scf
        MaxIter 1000
    end
    %output
        Print[ P_Hirshfeld] 1
    end
    %elprop
        Polar 1
    end
    %plots
        dim1 100
        dim2 100
        dim3 100
        Format Gaussian_Cube
        ElDens("{dens_file}");
        MO("{homo_file}", {homo_index}, 0);
        MO("{lumo_file}", {lumo_index}, 0);
    end
    %pal
        nprocs {nprocs}
    end
    *xyzfile 0 1 {xyz_filename}
    """

    functional = "B3LYP"
    dispersion_correction = True
    basis_set = "def2-tzvp"
    keepdens = True
    opt = True
    resp = False
    nprocs = 10
    
    def __init__(self):
        pass

    @classmethod
    def generate_input(cls, xyz_filepath, input_filepath, filename):

        homo_index, lumo_index = get_homo_lumo_from_xyz(xyz_filepath)
            
        keywords = f"!{cls.functional}"
        if cls.dispersion_correction:
            keywords += " D3BJ"
        keywords += (" " + cls.basis_set)
        if cls.opt:
            keywords += " Opt"
        if cls.resp:
            keywords += " chelpg"
            
        template_filled = cls.dft_template.format(
            keywords=keywords,
            dens_file=filename + ".dens.cube",
            homo_file=filename + ".homo.cube",
            lumo_file=filename + ".lumo.cube",
            homo_index=homo_index,
            lumo_index=lumo_index,
            nprocs=cls.nprocs,
            xyz_filename=filename + ".xyz"
        )
        
        with open(input_filepath, 'w') as file:
            file.writelines(template_filled)
        return template_filled

    @classmethod
    def set_opt(cls, opt):
        if not isinstance(opt, bool):
            print("Please provide True or False.")
            return()
        cls.opt = opt
        if opt == True:
            print("Geometry optimization will be executed.")
        if opt == False:
            print("No geometry optimization will be executed.")

    @classmethod
    def set_resp(cls, resp):
        if not isinstance(resp, bool):
            print("Please provide True or False.")
            return()
        cls.resp = resp
        if resp == True:
            print("Resp file will be generated.")
        if resp == False:
            print("No resp file will be generated.")

    @classmethod
    def set_functional(cls, functional):
        list_of_accepted_functionals = ["HF"]
        if functional in list_of_accepted_functionals:
            print(f"Functional set to {functional}")
            cls.functional = functional
        else:
            print("Functional not accepted. Using B3LYP as default.")

    @classmethod
    def set_dispersion_correction(cls, dispersion_correction):
        if not isinstance(dispersion_correction, bool):
            print("Please provide true or false")
        cls.dispersion_correction = dispersion_correction
        if dispersion_correction == True:
            print("Disperion correction will be applied.")

        if dispersion_correction == False:
            print("Dispersion correction will not be applied.")
            
    @classmethod
    def set_basis_set(cls, basis_set):
        list_of_accepted_basis_sets = ["6-31G*"]
        if basis_set in list_of_accepted_basis_sets:
            print(f"Basis set set to {basis_set}")
            cls.basis_set = basis_set
        else:
            print("Basis set not accepted. Using def2-tzvp as default.")

    @classmethod
    def set_nprocs(cls, nprocs):
        cls.nprocs = nprocs

    @classmethod
    def print_parameters(cls):
        print("Current parameters of DFT_input_generator:")
        print(f"Functional: {cls.functional}")
        print(f"Basis set: {cls.basis_set}")
        print(f"Geometry optimization: {cls.opt}")
        print(f"Dispersion correction: {cls.dispersion_correction}")