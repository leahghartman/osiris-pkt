import copy

"""
  ############################### How To USe This #######################
  TODO
"""

class NamelistContext(object):
  def __init__(self, name=None, label=None, desc=None, items=None, copy_template=None ):
    self.name = name
    self.label = label
    self.desc = desc
    self.items = items; self.items_dict={};
    self.copy_template = copy_template
    if self.copy_template is not None:
      self.items = copy.deepcopy( self.copy_template.items )
    for item in self.items:
      self.items_dict[item.name] = item
  def get( self, block_name=None, namelist_name=None):
    if block_name is not None and namelist_name is None:
      if block_name in self.items_dict:
        return self.items_dict[block_name]
      else:
         raise ValueError(" Error getting block '%s'. It is not found."
                             % ( block_name ))
    if block_name is None and namelist_name is not None:
      if namelist_name in  self.items_dict:
        return self.items_dict[namelist_name]
      else:
        raise ValueError(" Error getting namelist '%s'. It is not found."
                             % ( namelist_name ))
    if block_name is not None and namelist_name is not None:
      return self.items_dict[block_name].get( namelist_name )
    else:
      raise ValueError(" Error getting namelist.. both namelist_name and block_name were None")
  def replace( self, block_name=None, namelist_name=None, new_item=None):
    if block_name is not None:
      self.items_dict[block_name].replace( namelist_name, new_item)
    else:
      self.items[ self.items.index( self.items_dict[namelist_name] )] = new_item
      self.items_dict[namelist_name] = new_item
  def replace_all( self, namelist_name=None, new_item=None):
    if namelist_name in self.items_dict:
      self.replace(namelist_name=namelist_name, new_item=new_item )
    for item in self.items:
      if isinstance( item, NamelistOrderBlock):
        item.replace( namelist_name, new_item=new_item)
  def insert_after( self, block_name=None, after_namelist_name=None, new_item=None):
    if block_name is None:
      if after_namelist_name is not None and after_namelist_name in self.items_dict:
        self.items.insert(self.items.index( self.items_dict[after_namelist_name] ), new_item )
        self.items_dict[new_item.name] = new_item
      else:
        raise ValueError(" Error inserting after '%s'.'%s' is not found."
                             % ( after_namelist_name,self.name,after_namelist_name ))
    else:
      self.items_dict[block_name].insert_after( after_namelist_name, new_item )
class NamelistOrderEntry( object ):
  def __init__(self, name, required):
    self.name = name
    self.required  = required
class NamelistOrderBlock( object ):
  def __init__(self, name, repeats, how_to_tell_if_repeated, items):
    self.name = name
    self.items = items; self.items_dict={};
    self.repeats = repeats
    self.how_to_tell_if_repeated = how_to_tell_if_repeated
    self.items = items
    for item in items:
      self.items_dict[item.name] = item
  def get( self, namelist_name):
    if namelist_name in  self.items_dict:
      return self.items_dict[namelist_name]
    else:
      raise ValueError(" Error getting namelist '%s'. It is not found in block %s."
                           % ( namelist_name, self.name ))
  def replace( self, namelist_name, new_item ):
    if namelist_name in self.items_dict:
      self.items[ self.items.index( self.items_dict[namelist_name] )] = new_item
      self.items_dict[namelist_name] = new_item
  def insert_after( self, after_namelist_name, new_item):
    if after_namelist_name in self.items_dict:
      self.items.insert(self.items.index( self.items_dict[after_namelist_name] ), new_item )
    else:
      raise ValueError(" Error inserting after %s/%s. Block '%s' does not contain any item"
                           " with name '%s'" % ( self.name,after_namelist_name,self.name,after_namelist_name ))
"""
  Brief just-the-facts description of the ordering of the namelists... the order matters in Osiris.
    The structure of this list is:
      (namelist-name, context, is-required)

    namelist-name    : just the name given in the code to this namelist.
    namelist-context : one of the context defined above. It should be that the
                       combination of namelist-context+namelist_name is unique.
    is-required      : is this namelist required to always be present in a valid input file?
"""



"""
  Define the 'standard' ordering.. other simualtion contexts can use this as a basis.
"""
NL = NamelistOrderEntry
standard_template = NamelistContext( 
  name ='stndard_template', label="Standard Ordering Tempalte",
  desc ="Not meant to be used directly.. just copied", 
  items=[
      NamelistOrderBlock('General Simulation Parameters', repeats=False,  how_to_tell_if_repeated=None,
        items = [                 
          NL(                                       '/osiris/read_sim_options/nl_simulation',        True),
          NL(                                '/m_node_conf/read_nml_node_conf/nl_node_conf',         True),
          NL(                                               '/read_input_grid/nl_grid',              True),
          NL(                                '/m_time_step/read_nml_time_step/nl_time_step',         True),
          NL(                                    '/m_restart/read_nml_restart/nl_restart',           True),
          NL(                                        '/m_space/read_nml_space/nl_space',             True),
          NL(                                          '/m_time/read_nml_time/nl_time',              True),
        ]),
      NamelistOrderBlock('Electro-Magnetic Fields', repeats=False, how_to_tell_if_repeated=None, 
        items = [                         
          NL(                                            '/read_input_emf/nl_el_mag_fld',            False),
          NL(                            '/m_emf_bound/read_nml_emf_bound/nl_emf_bound',             True),
          NL(                              '/m_vdf_smooth/read_nml_smooth/nl_smooth',                False),
          NL(                            '/m_emf_diag/read_input_diag_emf/nl_diag_emf',              False)
        ]),
      NL(                                      '/read_input_particles/nl_particles',                 True),
      NamelistOrderBlock('Species', repeats=True, how_to_tell_if_repeated='particles.num_species', 
        items = [                         
          NL(                                        '/read_input_species/nl_species',               True),
          NL(                       '/m_species_udist/read_nml_spec_udist/nl_udist',                 False),
          NL(                              '/m_psource_std/read_input_std/nl_profile',               False),
          NL(                    '/m_psource_btfel/read_input_betatronfel/nl_betatron_fel',          False),
          NL(                     '/m_species_boundary/read_nml_spe_bound/nl_spe_bound',             True),
          NL(                                     '/m_piston/read_nml_pst/nl_piston',                False),
          NL(                                   '/read_input_diag_species/nl_diag_species',          False)
        ]),
      NamelistOrderBlock('Cathode', repeats=True, how_to_tell_if_repeated='particles.num_cathode', 
        items = [                         
          NL(                              '/m_cathode/read_input_cathode/nl_cathode',               True),
          NL(                                        '/read_input_species/nl_species',               True),
          NL(                       '/m_species_udist/read_nml_spec_udist/nl_udist',                 False),
          NL(                     '/m_species_boundary/read_nml_spe_bound/nl_spe_bound',             True),
          NL(                                   '/read_input_diag_species/nl_diag_species',          False)
        ]),
      NamelistOrderBlock('Neutrals', repeats=True, how_to_tell_if_repeated='particles.num_neutral', 
        items = [
          NL(                              '/m_neutral/read_input_neutral/nl_neutral',               True),
          NL(                                    '/m_cross/read_nml_cross/nl_num_ene',               False),
          NL(                                    '/m_cross/read_nml_cross/nl_cross',                 False),
          NL(                              '/m_psource_std/read_input_std/nl_profile',               False),
          NL(                    '/m_psource_btfel/read_input_betatronfel/nl_betatron_fel',          False),
          NL(                      '/m_diag_neutral/read_nml_diag_neutral/nl_diag_neutral',          False),
          NL(                                        '/read_input_species/nl_species',               True),
          NL(                     '/m_species_boundary/read_nml_spe_bound/nl_spe_bound',             True),
          NL(                                   '/read_input_diag_species/nl_diag_species',          False)
        ]),
      NamelistOrderBlock('Moving Neutral', repeats=True, how_to_tell_if_repeated='particles.num_moving', 
        items = [
          NL(                              '/m_neutral/read_input_neutral/nl_neutral_mov_ions',      True),
          NL(                              '/m_psource_std/read_input_std/nl_profile',               False),
          NL(                    '/m_psource_btfel/read_input_betatronfel/nl_betatron_fel',          False),
          NL(                      '/m_diag_neutral/read_nml_diag_neutral/nl_diag_neutral',          False),
          NL(                                        '/read_input_species/nl_species',               True),
          NL(                     '/m_species_boundary/read_nml_spe_bound/nl_spe_bound',             True),
          NL(                                        '/read_input_species/nl_species',               True),
          NL(                                   '/read_input_diag_species/nl_diag_species',          False),
          NL(                     '/m_species_boundary/read_nml_spe_bound/nl_spe_bound',             True),
          NL(                                   '/read_input_diag_species/nl_diag_species',          False)
        ]),
      NL(                            '/m_species_collisions/read_nml_coll/nl_collisions',            False),
      NamelistOrderBlock('Laser Pulses', repeats=True, how_to_tell_if_repeated='SELF', 
        items = [
          NL(                            '/m_zpulse_std/read_input_zpulse/nl_zpulse',                False),
          NL(              '/m_zpulse_mov_wall/read_input_zpulse_mov_wall/nl_zpulse_mov_wall',       False),
          NL(                    '/m_zpulse_point/read_input_zpulse_point/nl_zpulse_point',          False),
          NL(    '/m_zpulse_sliding_focus/read_input_zpulse_sliding_focus/nl_zpulse_sliding_focus',  False),
          NL(                      '/m_zpulse_wall/read_input_zpulse_wall/nl_zpulse_wall',           False),
          NL(                '/m_zpulse_speckle/read_input_zpulse_speckle/nl_zpulse_speckle',        False),
        ]),
      NamelistOrderBlock('Current', repeats=False, how_to_tell_if_repeated=None, 
        items = [
          NL(                                                           '/nl_current',               False),
          NL(                    '/m_current_diag/read_input_diag_current/nl_diag_current',          False)
        ]),
      NL(                    '/m_antenna_array/read_nml_antenna_array/nl_antenna_array',             False),
      NamelistOrderBlock('Antenna', repeats=True, how_to_tell_if_repeated="nl_antenna_array.num_antenna",
        items = [
          NL(                                '/m_antenna/read_nml_antenna/nl_antenna',               False)
        ])
    ])

"""
  ADD ANY NEW SIMULATION MODES INTO THIS LIST!
"""
namelist_ordering_contexts = {
  'std': NamelistContext( name='std', label="Standard", 
                   desc="The default Osiris Experience. When 'simulation.algorithm'"
                        " is not specified or then set to 'standard'",
                   copy_template=standard_template),
  'quasi-3D': NamelistContext( name='quasi-3D', label="Quasi-3D (Cylindirical Modes)", 
                   desc="When simulation.algorithm is set to 'quasi-3D'",
                   copy_template=standard_template),
  'fei': NamelistContext( name='fei', label="Fei Solver", 
                   desc="When simulation.algorithm is set to 'fei'",
                   copy_template=standard_template),
  'pgc': NamelistContext( name='pgc', label="Pondermotive Guiding Center", 
                   desc="When simulation.algorithm is set to 'pgc'",
                   copy_template=standard_template),
  'shear': NamelistContext( name='shear', label="Shear Mode", 
                   desc="When simulation.algorithm is set to 'shear'",
                   copy_template=standard_template),
  'qed':NamelistContext( name='qed', label="QED mode.", 
                   desc="When simulation.algorithm is set to 'qed'",
                   copy_template=standard_template),
  'radiat': NamelistContext( name='radiat', label="Radiation Mode", 
                   desc="When simulation.algorithm is set to 'radiat'",
                   copy_template=standard_template)
}

# NOW CUSTOMIZE the ordering of the particular simulation context!

#================ quasi-3D ===========================
namelist_ordering_contexts['quasi-3D'].replace(
    block_name = 'General Simulation Parameters',
    namelist_name = '/read_input_grid/nl_grid',
    new_item = NL( '/m_grid_cyl_modes/read_input_grid_cyl_modes/nl_grid', required = True) )
namelist_ordering_contexts['quasi-3D'].replace(
    block_name    = 'Electro-Magnetic Fields',
    namelist_name = '/read_input_emf/nl_el_mag_fld',
    new_item = NL( '/m_emf_cyl_modes/read_input_emf_cyl_modes/nl_el_mag_fld', required = False) )
namelist_ordering_contexts['quasi-3D'].replace_all(
    namelist_name = '/read_input_species/nl_species',
    new_item = NL( '/m_species_cyl_modes/read_input_species_cyl_modes/nl_species', required = False) )
namelist_ordering_contexts['quasi-3D'].insert_after(
    block_name = 'Current',
    after_namelist_name="/nl_current", 
    new_item = NL( '/m_current_cyl_modes/read_input_current_cyl_modes/nl_sp_filter', required = False) )


#================ Fei ===========================
# replaces 
#add /m_current_fei/read_input_fei/nl_sp_filter'
namelist_ordering_contexts['fei'].replace(
    block_name    = 'Electro-Magnetic Fields',
    namelist_name = '/read_input_emf/nl_el_mag_fld',
    new_item = NL(  '/m_emf_fei/read_input_fei/nl_el_mag_fld', required = True) )
namelist_ordering_contexts['fei'].insert_after(
    block_name = 'Current',
    after_namelist_name="/nl_current", 
    new_item = NL( '/m_current_cyl_modes/read_input_current_cyl_modes/nl_sp_filter', required = False) )

#================ PGC ===========================
namelist_ordering_contexts['pgc'].insert_after(
    block_name = 'Electro-Magnetic Fields',
    after_namelist_name="/m_emf_diag/read_input_diag_emf/nl_diag_emf", 
    new_item = NL( '/read_input_pgc/nl_pgc', required = False) )
namelist_ordering_contexts['pgc'].replace(
    block_name = 'Neutrals',
    namelist_name = '/m_neutral/read_input_neutral/nl_neutral',
    new_item = NL( '/m_neutral_pgc/read_input_neutral_pgc/nl_neutral', required = False) )                        
namelist_ordering_contexts['pgc'].replace(
    block_name = 'Moving Neutral',
    namelist_name = '/m_neutral/read_input_neutral/nl_neutral_mov_ions',
    new_item = NL( '/m_neutral_pgc/read_input_neutral_pgc/nl_neutral_mov_ions', required = False) )                        

#================ Shear ===========================
namelist_ordering_contexts['shear'].insert_after(
    block_name = 'Electro-Magnetic Fields',
    after_namelist_name="/m_emf_diag/read_input_diag_emf/nl_diag_emf", 
    new_item = NL( '/m_emf_shear/read_input_shear/nl_shear', required = False) )

#================ Radio ===========================
namelist_ordering_contexts['radiat'].replace(
    namelist_name = '/read_input_particles/nl_particles',
    new_item = NL( '/m_particles_rad/read_input_rad/nl_particles', required = True) )
namelist_ordering_contexts['radiat'].insert_after(
    block_name = 'Species',
    after_namelist_name="/read_input_diag_species/nl_diag_species", 
    new_item = NL( '/class_CartDet/read_input_detector/nl_cartesian', required = False) )
namelist_ordering_contexts['radiat'].insert_after(
    block_name = 'Species',
    after_namelist_name="/read_input_diag_species/nl_diag_species", 
    new_item = NL( '/class_SpheDet/read_input_detector/nl_spherical', required = False) )


#================ Qed ===========================
"""
# replaces /m_particles_qed/read_input_qed/nl_particles
  adds a species group called a QED group. it consists of
   - species block 'electrons' (stanrard species with megere info at the end)
   - species block 'psoitrons'  (stanrard species with megere info at the end)
   - new kinf of species block 'photons' with are comprised of:
        NL(                    '/m_photons_class/read_input_qed_photons/nl_photons',                   False),
        NL(                       '/m_species_udist/read_nml_spec_udist/nl_udist',                     False),
        NL(                       '/m_photons_class/read_nml_phot_bound/nl_phot_bound',                False),
        NL(                        '/m_photons_class/read_nml_phot_diag/nl_diag_photons',              False),
        NL(                          '/m_mergedel/read_input_merge_info/nl_merge_info',                True),
        NL(               '/m_particles_qed/read_input_qed_collide_info/nl_qed_collide_info',          True),
  
"""
"""
# QED gets rid of cathods, Neutrals.
namelist_ordering_contexts['qed'].remove(
    namelist_name = "EMF")
namelist_ordering_contexts['qed'].remove(
    namelist_name = "Cathode")
namelist_ordering_contexts['qed'].remove(
    namelist_name = "Neutrals")
namelist_ordering_contexts['qed'].remove(
    namelist_name = "Moving Neutral")
# namelist_ordering_contexts['qed'].replace(
#     namelist_name = '/read_input_particles/nl_particles',
#     new_item = NL( '/m_particles_qed/read_input_qed/nl_particles', required = True) )
# after the particles
#  make a special species block for electrons
electron_block = copy.deepcopy( namelist_ordering_contexts['qed'].get(block_name="Species") )
positron_block = copy.deepcopy( namelist_ordering_contexts['qed'].get(block_name="Species") )
"""
