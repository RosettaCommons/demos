# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
import os,sets,sys,re,gc

sys.setrecursionlimit(10000)

from pygccxml import declarations
import pygccxml
import pyplusplus
from pyplusplus.module_builder import call_policies

sys.path.append("../../")
from pyplusplus_utils import *


def make_mb(srcfiles):
  if type(srcfiles) is type(""): srcfiles = [ srcfiles ]
  path="/Users/sheffler/svn/branches/mini"
  mb = pyplusplus.module_builder.module_builder_t ( srcfiles
        , gccxml_path="/usr/local/bin"
        , working_directory = path+"/demo/python/by_hand/core"
        , include_paths   = [ path+'/src', path+'/src/platform/linux' ]
        , indexing_suite_version=1
        , define_symbols=[]
        , cache =  pygccxml.parser.file_cache_t( "/tmp/core_all.gccxml.cache" )
        )
  return mb


del mb
gc.collect()
mb = make_mb("include/mini_all.hh")


def isinc(x):
  return not x._ignore

def is_refcount(c):
  return "utility::pointer::ReferenceCount [class]" in [str(b) for b in get_class_bases(c)]

def postproc(code):
  code = re.sub(r'<unsigned([,>])',r'<std::size_t\1',code)
  return code

def excludes(mb):
  mb.class_("ContextGraph").exclude()
  mb.class_("DistanceNode").exclude()
  mb.class_("DistanceEdge").exclude()
  mb.class_("DistanceGraph").exclude()
  mb.class_("PointGraphVertexData").var("NUM_EDGES_TO_RESERVE").exclude()
  mb.class_("FoldTree").var("PEPTIDE").exclude()
  mb.class_("FoldTree").mem_fun("count_fixed_residues").exclude()
  mb.class_("ScoreFunction").mem_fun("setup_for_scoring").exclude()
  mb.class_("EnergyGraph").mem_fun("energy_exists").exclude()
  mb.free_funs("fullatom_residue_pair_energy").exclude()
  mb.free_funs("score_function_hb_eval_type_weight").exclude()
  mb.free_funs("hbond_energy").exclude()
  mb.namespace("kinematics").mem_funs("get").exclude()


def ischildof(a,b):
  while a and a is not b: a = a.parent
  return a is b


#
#
#   FIX held types!!! ContextIndependentOneBodyEnergy_wrapper
#
#

""" must add this to c'tors for energy method wrappers!!!!
#using namespace core::scoring; APL: Will, how do I fix this?
#add_score_type(python);
"""

allprotocol = "RigidBodyPerturbMover ProtocolOptions RottrialOptions RigidBodyMover PackMover ShearMover SmallMover MoveSet Mover MonteCarlo ClassicRelax DockingHighRes Score ScoreManager ProtocolOptions RottrialOptions"

corecls = 'Pose ResidueType AtomType ResidueTypeSet AtomTypeSet EMapVector DOF_ID AtomTreeMinimizer PackerTask Residue Ramachandran MoveMap ScoreFunction PackerTask MinimizerOptions EnergyMethod ContextIndependentOneBodyEnergy ContextIndependentTwoBodyEnergy '

def addcode(mb):
  while mb.registrations_code_head: mb.registrations_code_tail.pop()
  while mb.registrations_code_tail: mb.registrations_code_tail.pop()
  while mb.declarations_code_head: mb.declarations_code_tail.pop()
  while mb.declarations_code_tail: mb.declarations_code_tail.pop()
  mb.add_declaration_code("""
  #include "boost/python/suite/indexing/indexing_suite.hpp"
  #include "boost/python/suite/indexing/vector_indexing_suite.hpp"
  #include "boost/python/suite/indexing/map_indexing_suite.hpp"\n""")
  mb.add_declaration_code("""void init_just_db_for_python_tests(std::string const & dbloc) {
    int argc = 3;
    char ** argv = new char*[4];
    argv[0] = new char[999];
    argv[1] = new char[999];
    argv[2] = new char[999];
    argv[3] = new char[999];
    std::strcpy(argv[0],"dummy");
    std::strcpy(argv[1],"-database");
    std::strcpy(argv[2],dbloc.c_str());
    std::strcpy(argv[3],"-constant_seed");
    core::init( argc, argv );
  }
  """)
  mb.add_registration_code('bp::def("init",&init_just_db_for_python_tests)')
  mb.add_registration_code("""using namespace utility::pointer;
  using namespace protocols::moves;
  bp::implicitly_convertible< owning_ptr<SmallMover>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<ShearMover>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<RigidBodyPerturbMover>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<RigidBodyMover>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<PackMover>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<MoveSet>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<MoveSet_wrapper>, owning_ptr<MoveSet> >();
  bp::implicitly_convertible< owning_ptr<Mover_wrapper>, owning_ptr<Mover> >();
  bp::implicitly_convertible< owning_ptr<PackerTask_wrapper>, owning_ptr<core::pack::task::PackerTask> >();
  bp::implicitly_convertible< owning_ptr<ScoreFunction_wrapper>, owning_ptr<core::scoring::ScoreFunction> >();
  bp::implicitly_convertible< owning_ptr<ContextIndependentOneBodyEnergy_wrapper>, owning_ptr<core::scoring::methods::ContextIndependentOneBodyEnergy> >();
  bp::implicitly_convertible< owning_ptr<core::scoring::methods::ContextIndependentOneBodyEnergy>, owning_ptr<core::scoring::methods::EnergyMethod> >();
  bp::class_<     std::vector <std::string> > ( "VecStr").def(bp::vector_indexing_suite< std::vector<std::string> >());
  bp::class_< utility::vector1<std::string> > ("Vec1Str").def(bp::vector_indexing_suite< utility::vector1<std::string> >());
  bp::class_<     std::vector <int> > ("Vec1Int").def(bp::vector_indexing_suite< std::vector<int> >() );
  bp::class_< utility::vector1<int> > ( "VecInt").def(bp::vector_indexing_suite< utility::vector1<int> >() );
  bp::class_<     std::vector <bool> > ("Vec1bool").def(bp::vector_indexing_suite<      std::vector<bool> >() );
  bp::class_< utility::vector1<bool> > ( "Vecbool").def(bp::vector_indexing_suite< utility::vector1<bool> >() );
  """)

addcode(mb)


def wrap_relax(pcls,ccls,ffuns):
  core = mb.namespace("core")
  prtl = mb.namespace("protocols")
  mb.global_ns.exclude()
  for f in mb.calldefs().declarations:
    f.create_with_signature = True
    if f.call_policies is None:
      f.call_policies = call_policies.return_value_policy( call_policies.return_by_value )
  print 'line 148'
  for f in ffuns.split(): mb.free_funs(f).include()
  for cname in pcls.split():
    print 'including',cname
    prtl.class_(cname).include()
  for cname in ccls.split():
    print 'including',cname
    core.class_(cname).include()
  core.namespace("conformation").class_("Atom").include()
  for ns in (core,prtl):
    ns.enums().include()
    for c in ns.classes():
      # if c.name == 'ScoreFunction':
      #   print '\n\n\n\n\n\n',is_refcount(c), c.is_abstract, len(c.is_wrapper_needed()),'\n\n\n\n\n\n'
      # if is_refcount(c):
      #   if c.is_abstract or len(c.is_wrapper_needed()) > 0:
      #     c.held_type = "utility::pointer::owning_ptr< "+c.wrapper_alias+" >"
      #   else: c.held_type = "utility::pointer::owning_ptr< "+c.decl_string+" >"
      for p in c.private_members: p.exclude()
      for p in c.protected_members: p.include()
  #
  print "168"
  for o in mb.global_ns.operators("=",allow_empty=True).declarations: o.exclude()
  mb.class_("PackerTask").mem_fun("clone").include()
  mb.class_("EnergyMethod").mem_fun("clone").include()
  mb.class_("DockingHighRes").mem_funs("call_pack").exclude()
  for cname in "MoveMap Residue Pose".split():
    while core.class_(cname).registration_code: core.class_(cname).registration_code.pop()
  core.class_('MoveMap').add_code('def("__iter__", bp::range( &core::kinematics::MoveMap::dof_begin, &core::kinematics::MoveMap::dof_end ))')
  core.class_("Residue").add_code('def("__iter__", bp::range( &core::conformation::Residue::atom_begin, &core::conformation::Residue::atom_end))')
  core.class_("Pose"   ).add_code('def("__iter__", bp::range( &core::pose::Pose::res_begin, &core::pose::Pose::res_end))')
  #
  core.class_("ScoreFunction").mem_fun("setup_for_scoring").exclude()
  #core.class_("EnergyMethod").mem_fun("add_score_type").include() -- APL: Will, how do I fix this?
  #
  for f in mb.calldefs().declarations:
    f.create_with_signature = True
    if f.call_policies is None:
      if not hasattr(f,'return_type') or not f.return_type: continue
      if f.return_type.decl_string.count("const"):
        f.call_policies = call_policies.return_value_policy( call_policies.copy_const_reference )
      else:
        f.call_policies = call_policies.return_value_policy( call_policies.copy_non_const_reference )
  #
  #
  name = "test_relax"
  mb.build_code_creator( '_'+name )
  #
  for ns in (core,prtl):
    for c in ns.classes():
      if is_refcount(c):
        if c.is_abstract:# or len(c.is_wrapper_needed()) > 0:
          c.held_type = "utility::pointer::owning_ptr< "+c.wrapper_alias+" >"
        else: c.held_type = "utility::pointer::owning_ptr< "+c.decl_string+" >"
  #mb.class_("PackerTask").held_type = "utility::pointer::owning_ptr<core::pack::task::PackerTask>"
  #
  print len(filter(isinc,mb.decls()))
  creators = pyplusplus.code_creators.make_flatten(mb.code_creator.creators)
  ccc =       filter( lambda x: isinstance(x,pyplusplus.code_creators.mem_fun_pv_t), creators )
  ccc.extend( filter( lambda x: isinstance(x,pyplusplus.code_creators.mem_fun_pv_wrapper_t), creators ) )
  ccc.extend( filter( lambda x: isinstance(x,pyplusplus.code_creators.mem_fun_v_t), creators ) )
  ccc.extend( filter( lambda x: isinstance(x,pyplusplus.code_creators.mem_fun_v_wrapper_t), creators ) )
  ccc = filter( lambda x: x.parent.declaration.name in ("ClassicRelax") ,ccc)
  print [c.parent.creators.remove(c) for c in ccc if c.declaration.name == "clone"]
  for c in "ClassicRelax".split():
    mb.class_(c).mem_fun("clone").exclude()
  # for cc in filter( lambda x: isinstance(x,pyplusplus.code_creators.class_wrapper_t), creators):
  #   if is_refcount(cc.declaration):
  #      cc.class_creator._set_held_type("utility::pointer::owning_ptr< "+cc.wrapper_alias+" >")
  # ns = namespace_from_full_path(mb,'core')
  # not_in_ns = filter( lambda x: not ischildof(x.declaration,ns) , ccc )
  # for x in not_in_ns: x.parent.creators.remove(x)
  # print len(not_in_ns)
  mb.write_module( './wrap/wrap_'+name+'.cc' )
  f = open('./wrap/wrap_'+name+'.cc')
  s = postproc(f.read())
  f.close()
  f = open('./wrap/wrap_'+name+'.cc','w')
  f.write(s)
  f.close()
  return len(filter(isinc,mb.decls()))


wrap_relax(allprotocol,corecls,"pose_from_pdb")

















