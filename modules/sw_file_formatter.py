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
        cls.opt = opt
        if opt == True:
            print("Geometry optimization will be executed.")
        if opt == False:
            print("No geometry optimization will be executed.")

    @classmethod
    def set_functional(cls, functional):
        list_of_accepted_functionals = []
        if functional in list_of_accepted_functionals:
            cls.functional = functional
        else:
            print("Functional not accepted. Using B3LYP as default.")

    @classmethod
    def set_nprocs(cls, nprocs):
        cls.nprocs = nprocs