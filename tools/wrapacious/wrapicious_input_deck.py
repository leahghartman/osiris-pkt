import glob

# Specific configuration for Osiris...
#   Let's add all files that have benn preprocessed as part of the compilation process of Osiris itself.
SourceFilesGlobs = [ "./preprocessed/*.f90", "./preprocessed/*.f03"]

# Plus exclude a few things we won't need in the wrapper.
SourceFilesExcludes = [ '*/ioperf.f90', 
                        '*-meta*', 
                        '*os-psource-def*',
                        '*os-spec-diag-def*',  
                        '*os-spline.f90*',
                        '*zdf-interface*',
			                  '*memory*',
                        '*stringutil*',
                        '*fparser*']
SourceFilesList = []


def WrapaciousConfiguration():
  return {

    # Various directories...
    "out_dir":        "./processed_wrapper",
    "bin_dir":        "./",
    "in_dir":         "./",
    "out_dir_python": "./processed_wrapper/python",
    "temp_dir":       "./processed_wrapper/cache",

    # If True, always parse all input source code file from scratch - ignoring any cached results from previous runs (like doing 'make clean' every time). 
    #   Not recommended for production. Great for debugging.
    "always_read_input": False,

    # Used when debugging so that the wrapper is always generated (and the wrapping code always executed), 
    #   even if the input source code (i.e. stuff we are wrapping) is unchanged.
    "always_write_output": False,                        

    # you prolly want this to be True unless you have specific reasons otherwise.
  	"removePreprocessorDirectives": True,

    #
    # Python Interface Options
    #
    "fortran_string_class_name": "_str",
    
    # Use our previously specified source files and officially add them into the options.
    "SourceFilesList":     SourceFilesList,
    "SourceFilesGlobs":    SourceFilesGlobs,
    "SourceFilesExcludes": SourceFilesExcludes,
    
    "final_python_filename": "pyosiris",
    "final_wrapper_sharedlib_filename": "pyosiris_lib",
    "alterationList": [
        { "type": "use", "file": "os-emf.f90", "global":True, "module": None, "routine":"init_emf", "match":"vdf_memory", "delete":True, "replace":None },
        { "type": "derived type var", "file": "os-spec-define.f90", "module": None, "derived type": "t_species", "match":"events", "delete":True, "replace":None },
    ],
    "array options": {
        "index order": "c",                     # valid settings are 'fortran' or 'c'. The default is 'c'
        "base one default": "False",            # If the arrays, when returned, default to index with a base of 1
    },
    "ignore items": [   'fortran_string___new',
                        'getenv_f',
                        'gethostname_f',
                        'mkdir_f',
                        'chdir_f',
                        'getcwd_f',
                        'strerror_f',
                        'symlink_f',
                        'remove_f',
						'alloc_1d_str',
						'alloc_2d_str',
						'alloc_3d_str',
						'alloc_4d_str',
						'alloc_bound_1d_str',
						'alloc_bound_2d_str',
						'alloc_bound_3d_str',
						'alloc_bound_4d_str',
						'freemem_1d_str',
						'freemem_2d_str',
						'freemem_3d_str',
						'freemem_4d_str',
						'li_ene',
                        'li_cross',
                        'get_str',          # good example oof char buffer versus string.. ;dont handle it correctively  yet
            'avail_report_quants_current_cyl_modes',  # Name too long
            'init_report_quants_current_cyl_modes',    # Name Too Long
            'update_boundary_bfld_emf_bound_cyl_modes',
            'update_boundary_bfld_emf_bound_cyl_modes_2',
            'update_boundary_efld_emf_bound_cyl_modes',
            'update_boundary_efld_emf_bound_cyl_modes_2',
            'operator(==)',
            'operator(/=)',
            'assignment(=)',
            'm_fft'
    ],
  }

  
