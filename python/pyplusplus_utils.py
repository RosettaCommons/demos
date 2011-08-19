# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
import os,sets,sys,re

from pygccxml import declarations
import pygccxml

import pyplusplus
from pyplusplus import code_creators, module_builder
from pyplusplus.module_builder import call_policies

ALIGN_OUTPUT = False
sys.path.append(os.environ["HOME"]+"/scripts")
try:
  from aligntext import alignMultiple
  ALIGN_OUTPUT = True
except:
  pass

from math import ceil
sys.path.append("~/scripts")



############################################################################
####################### utilities ##########################################

class CodeCreatorHolder(object):
  """docstring for CodeCreatorHolder"""
  def __init__(self, mb, code_creator=None):
    self.mb = mb
    if code_creator:
      self.code_creator = code_creator
    else:
      self.code_creator = mb.code_creator
    self.creator_type_names = filter(lambda x: x.endswith("_t"),dir(pyplusplus.code_creators))

  def creator_types_by_regex(self,regex,reopts=False,pad=True):
    if pad:
      regex = ".*%s.*"%regex
    type_names = filter( re.compile(regex,reopts).match, self.creator_type_names )
    return tuple([getattr(pyplusplus.code_creators,x) for x in type_names])

  def find(self,name=None,code=None,decider_func=None,reopts=re.MULTILINE|re.DOTALL,decl=None \
          ,klass=None,types=pyplusplus.code_creators.declaration_based_t):
    self.creators = pyplusplus.code_creators.make_flatten( self.code_creator.creators )
    if type(types) is type(''):
      types = self.creator_types_by_regex(types,reopts)
    creators = self.creators
    if types:
      creators = filter(lambda x: isinstance(x,types),creators)
    if decider_func:
      creators = filter(decider_func,creators)
    if code: # if string make re, if re use it
      f = lambda x: re.compile(code,reopts).match(x.create())
      creators = filter(f,creators)
    if name: # if string make re, if re use it
      f = lambda x: re.compile(name,reopts).match(x.declaration.name)
      creators = filter(f,creators)
    if decl: # if string make re, if re use it
      f = lambda x: re.compile(decl,reopts).match(x.declaration.decl_string)
      creators = filter(f,creators)
    if klass: # if string make re, if re use it
      f = lambda x: re.compile(decl,reopts).match(x.declaration.parent.decl_string)
      creators = filter(f,creators)
    return creators


def parse_full_name(full_name):
  pth = []#filter(len,full_class_name.split("::"))
  regex = re.compile("^([^<>, :]*)::")
  while regex.match(full_name):
    ns = regex.search(full_name).group(1)    # print ns
    if ns:
      pth.append( ns )
    full_name = regex.sub("",full_name)    # print full_name
  pth.append(full_name)
  return pth

def class_from_full_path(mb,fullname):
  path = parse_full_name(fullname)
  nspath,cname = path[:-1],path[-1]
  ns = mb.global_ns
  for nsname in nspath:
    for tmp in ns.namespaces(nsname):
      if tmp.parent is ns:
        ns = tmp
  return ns.class_(cname)

def namespace_from_full_path(mb,fullname):
  nspath = parse_full_name(fullname)
  ns = mb.global_ns
  for nsname in nspath:
    for tmp in ns.namespaces(nsname):
      if tmp.parent is ns:
        ns = tmp
  return ns

def print_class(c,out=sys.stdout):
  """docstring for print_class"""
  print >>out,c
  for pm in c.public_members:
    print >>out,'   ',pm

def rmns(s):
  return re.sub("(\w+::)+","",s)

def print_classes(mb,cnames,out=None):
  """docstring for print_classes"""
  if out is None:
    out = open('class_desc.txt','w')
  pmre = re.compile("""( .+? (?:<.*>)? ) [ ] ( (?:const[ ])? )  ( (?:[&*][ ])? )  (\S+)  [(] (.*?) [)]  ( (?:[ ]const)? )  [ ] \[.*?\] $ \n# brackets """,re.VERBOSE)
  for cname in cnames:
    p = parse_full_name(cname)
    basens = '::'.join(p[:-1])
    basecname = p[-1]
    c = class_from_full_path(mb,cname)
    s = "/"*(40-int(ceil(len(cname)/2)))+' '+cname+' '+"/"*(40-(len(cname)/2))
    print >>out, '/'*len(s)
    print >>out, s
    print >>out, '/'*len(s)
    print >>out, "namespace "+basens+"{"
    print >>out, "class", basecname, ": ( inherits?? )\n{"
    ret,rettype,retconst,retref,name,args,const = [],[],[],[],[],[],[]
    funstr = ""
    for pm in c.public_members:
      s = rmns(str(pm))
      if   s.endswith("[typedef]"):
        continue
      elif s.endswith("[constructor]"):
        continue
      elif s.endswith("[copy constructor]"):
        continue
      elif s.endswith("[class]"):
        continue
      elif s.endswith("[variable]"):
        continue
      elif s.endswith("[member function]"):
        a,b,c,d,e,f= pmre.findall( s )[0]
        t0re = "(?:%(template)s)"
        t1re = "(?:%(template)s(?:<"+t0re%dict(template="\w+\**&?")+"(?:, "+t0re%dict(template="\w+\**&?")+")*[ ]?>)?)"
        t2re = "(?:%(template)s(?:<"+t1re%dict(template="\w+\**&?")+"(?:, "+t1re%dict(template="\w+\**&?")+")*[ ]?>)?|"+t1re+")"
        t3re = "(?:%(template)s(?:<"+t2re%dict(template="\w+\**&?")+"(?:, "+t2re%dict(template="\w+\**&?")+")*[ ]?>)?|"+t2re+")"
        a = re.sub(t3re%dict(template="allocator"),'',a)
        a = re.sub(",\s+",",",a)
        a = a.replace("__normal_iterator","iter")
        a = a.replace("_iterator","_iter")
        a = a.replace("_const_iterator","_citer")
        ret.append(a.strip()+' '+b.strip()+' '+c.strip())
        rettype.append(a.strip())
        retconst.append(b.strip().ljust(5))
        retref.append(c.strip().ljust(1))
        name.append(d.strip())
        args.append(e.strip())
        const.append(f.strip().ljust(5))
    if name:
      m        = max( [len(s.strip()) for s in ret ])
      ret      = [ s.strip().rjust(m) for s in ret ]
      m        = max( [len(s.strip()) for s in rettype ])
      rettype  = [ s.strip().ljust(m) for s in rettype ]
      m        = max( [len(s.strip()) for s in retconst ])
      retconst = [ s.strip().ljust(m) for s in retconst ]
      m        = max( [len(s.strip()) for s in retref ])
      retref   = [ s.strip().ljust(m) for s in retref ]
      m        = max( [len(s.strip()) for s in name ])
      name     = [ s.strip().rjust(m) for s in name ]
      m        = max( [len(s.strip()) for s in const ])
      const    = [ s.strip().ljust(m) for s in const ]
      argl = max( [ max([len(s.strip()) for s in a.split(',') ]) for a in args ] )
      # print args
      # print argl
      for i in range(len(name)):
        head = '  '+ret[i] #' '.join((rettype[i],retconst[i],retref[i]))
        a = args[i].split(',')
        if len(a) > 1:
          l = len(head) + len(name[i]) + 2
          a = ( "\n" + " "*l + ',' ).join([s.ljust(argl) for s in a])
          funstr += ( head + ' ' + name[i]+' ( ' + a        +"\n"+" "*l  + ') '+const[i] ).rstrip() + ";\n"
        elif a[0] == "":
          funstr += ( head + ' ' + name[i]+' ( ' + ' '*argl   + ' ) ' + const[i] ).rstrip() + ";\n"
        else: # one arg
          funstr += ( head + ' ' + name[i]+' ( ' + a[0].strip().ljust(argl) +' ) ' + const[i] ).rstrip() + ";\n"
    if 0:
      funlines = funstr.split("\n")
      l = 999
      for line in funlines:
        if not line:
          continue
        gpst = line.find("  ")
        i = gpst
        while line[i] == ' ': i += 1
        nsp = i-gpst+2
        if nsp < l: l = nsp
      # for i,line in enumerate(funlines):
      #   if not line:
      #     continue
      #   gpst = line.find(" "*l)
      #   print l
      #   print line
      #   print line[:gpst],'.......',line[gpst+l:]
      #   funlines[i] = line[:gpst]+line[gpst+l:]
      funstr = "\n".join(funlines)
    print >>out, funstr+"}"
    print >>out,"}"
    print >>out

def get_classes_from_ns_path(full_class_name,starting_ns):
  pth = parse_full_name( full_class_name )
  print pth
  namespace = starting_ns
  for ns in pth[:-1]:
    print "INCLUDE CLASSES: looking for namespace",ns
    namespace = namespace.namespace(ns,recursive=False)
  classes = []
  ns_or_class = pth[-1]
  # if full_class_name path ends in namespace, add all classes in it
  for ns in namespace.namespaces(ns_or_class,allow_empty=True).declarations:
    classes.extend(ns.classes())
  # if full_class_name ends in class(s), find all in preceeding NS and add
  for c in namespace.classes(ns_or_class,allow_empty=True).declarations:
    classes.append(c)
  return classes


def get_classes_and_bases_from_ns_path(full_class_name,starting_ns):
  classes_and_bases = []
  for c in get_classes_from_ns_path(full_class_name,starting_ns):
    classes_and_bases.extend( get_class_bases(c) )
  return classes_and_bases


def get_type_base(arg):
  b = arg
  while hasattr(b,'base'): b = b.base
  assert hasattr(b,'declaration')
  return b.declaration


def get_class_bases(argcls,baseslist=None):
  if hasattr(argcls,'bases'):
    for bi in argcls.bases:
      assert hasattr(bi,"related_class")
      baseslist = get_class_bases(bi.related_class,baseslist)
  if baseslist is None:
    baseslist = []
  baseslist.append(argcls)
  return baseslist


def get_closest_superclass_returner(decls):
  bases = {}
  for d in decls:
    bases[d] = get_class_bases(get_type_base(d.declaration.return_type))
  blen = [len(b) for b in bases.values()]
  if sum([l == max(blen) for l in blen]) == 1:
    for d in decls:
      if len(bases[d]) == max(blen):
        return d,bases[d][-1]
  print "WARNING: couldn't decide which of ambiguous member functions to keep!"
  return None


def prune_ambiguous_function_overloads(class_creator,dups,conservative=True):
  keeper,super_class = get_closest_superclass_returner(dups)

  print dups[0].declaration.name
  impl_class  = class_creator.declaration.decl_string
  return_type = super_class.decl_string
  print return_type
  print impl_class
  if conservative and return_type != impl_class:
    print "WARNING: ignoring our choice of ambiguous function:"
    print "        ",keeper.declaration
    print "         because its return_type base isnt the class type:"
    print "        ",impl_class,"that were looking for"
    print "         this is a HACK... doesnt handle Keys and such"
    print "         being conservative about what to wrap so things will compile"
    print "         also, get_closest_superclass_returner doesnt seem"
    print "         to work on keys"
    keeper = None
  else:
    print "      keeping function which returns",return_type
  for d in dups:
    if d != keeper:
      class_creator.creators.remove(d)

  if len(class_creator.creators) == 0:
    print "WARING, class has no creators left!"

  # keepers = []
  # code_create_class = class_creator.declaration.decl_string
  # names = [d.declaration.parent.decl_string for d in dups]
  # exact_match = False
  # if code_create_class in names:
  #   exact_match = True
  # for d in dups:
  #   method_defined_class = d.declaration.parent.decl_string
  #   if not exact_match:
  #     method_defined_class += '_'
  #   if method_defined_class == code_create_class:
  #     keepers.append(d)
  # if len(keepers) == 1:
  #   keeper = keepers[0]
  #   for d in dups:
  #     if d is not keeper:
  #       d.parent.creators.remove(d)
  # else:
  #   print "WARNING: could not disambiguate "
  #   print code_create_class
  #   for k in dups:
  #     print k.declaration
  #   if len(keepers):
  #     print "I liked these:"
  #     for k in keepers:
  #       print k.declaration.parent.decl_string
  #
  #   sys.exit()


def fund_ambiguous_function_overloads(code_creator):
  duplicates = {}
  for cc in code_creator:
    t = pyplusplus.code_creators.calldef
    if not isinstance( cc, member_function_code_creators() ):
      continue
    d = cc.declaration
    k = (cc.parent,d.name,tuple(d.argument_types),d.has_const) # creator class, name, args, const
    if not k in duplicates:
      duplicates[k] = []
    duplicates[k].append(cc)
  for k in duplicates.keys():
    if len(duplicates[k]) == 1:
      del duplicates[k]
  return duplicates


def print_namespace_tree(mb,ns,depth=0,spacer="  "):
  print depth*spacer,ns.name
  for c in ns.namespaces(allow_empty=True):
    print_namespace_tree(mb,c,depth+1)


def make_alias(s):
  r = s
  regex = re.compile("""([,<> ]|^)
                            (?: [:][:])?
                            (?: [^:,<>.]+
                             [:][:] )+
                             ([^:,<>]+)""",re.VERBOSE)
  r = regex.sub("\g<1>\g<2>",s,re.VERBOSE)

  r = re.sub("[<>]","_",r)
  r = re.sub("[ ,]+","_",r)
  r = re.sub("[&*]","",r)
  while len(r) and r[-1] in "-_":
    r = r[:-1]
  return r



def basetype(typestr):
  """docstring for basetype"""
  typestr = typestr.split("::")[-1].strip()
  return re.sub("(?:const|[&]|<.*?>)","",typestr).strip()


def make_function(fun):
  ret = get_fun_return_type(fun)
  fun_name = declarations.full_name(fun)
  fun_sig  = fun.decl_string
  fstr = "( %(fun_sig)s )( & %(fun_name)s )"%vars()
  if '&' in ret or '*' in ret:    ###!!!!!!! this is where call policies are set in this script
    if 'const' in ret:
      cp = "bp::return_value_policy<bp::copy_const_reference>"
    else:
      cp = "bp::return_value_policy<bp::copy_non_const_reference>"
    return """bp::make_function( %(fstr)s , %(cp)s() )""".strip()%vars()
  else:
    return fstr



############################################################################
#################### libRosetta property recognizer ########################
class RosettaPropertyRecognizer(pyplusplus.decl_wrappers.properties.name_based_recognizer_t):
  def __init__(self):
    pyplusplus.decl_wrappers.properties.property_recognizer_i.__init__(self)

  def prefixes( self ):
    return [ ('','') ]

  def create_property( self, fget, fset ):
    if not fget.name == fset.name:
      return
    print "READ WRITE PROP:",fset.name
    # print "         getter:",SimpleFun(fget)
    # print "         setter:",SimpleFun(fset)
    prop_t = pyplusplus.decl_wrappers.properties.property_t
    return prop_t( fget.name, fget, fset )

  def is_accessor( self, mem_fun ):
      if mem_fun.ignore:
          return False
      if mem_fun.access_type != 'public':
          return False
      if mem_fun.has_static:
          return False #TODO: should be supported
      if mem_fun.virtuality == declarations.VIRTUALITY_TYPES.PURE_VIRTUAL:
          return False
      return True

  def is_getter( self, mem_fun ):
      if mem_fun.arguments:
          return False
      if declarations.is_void( mem_fun.return_type ):
          return False
      if not mem_fun.has_const:
          return False
      # print "think is getter:", SimpleFun(mem_fun)
      # print "                ", mem_fun.parent.name
      for o in mem_fun.overloads:
        if len(o.arguments)==0 and not declarations.is_void(o.return_type) and o.has_const:
          # print "overloaded also looks like getter:", SimpleFun(o)
          return False
      return True

  def is_setter( self, mem_fun ):
      if len( mem_fun.arguments ) != 1:
          return False
      if not declarations.is_void( mem_fun.return_type ):
          return False
      if mem_fun.has_const:
          return False
      if mem_fun.overloads:
        for o in mem_fun.overloads: # if any overloads also look like setters, don't bother
          if len( o.arguments ) == 1 and declarations.is_void( o.return_type ) and not o.has_const:
            return False
      # print "think is setter:", SimpleFun(mem_fun)
      # print "                ", mem_fun.parent.name
      return True

  def __get_accessors( self, mem_funs ):
      getters = []
      setters = []
      for mem_fun in mem_funs:
          if not self.is_accessor( mem_fun ):
              continue
          elif self.is_getter( mem_fun ):
              getters.append( mem_fun )
          elif self.is_setter( mem_fun ):
              print "FOUND SETTER!",mem_fun.name,mem_fun.parent.name
              setters.append( mem_fun )
          else:
              continue
      return ( getters, setters )



class RosettaPropertyFinder(pyplusplus.decl_wrappers.properties.properties_finder_t):
  """docstring for RosettaPropertyFinder"""
  def __init__(self, cls, recognizer=None, exclude_accessors=True ):
    self.super = pyplusplus.decl_wrappers.properties.properties_finder_t
    self.super.__init__(self,cls,recognizer,exclude_accessors)

  def find_properties( self, getters, setters, used_getters, used_setters ):
    properties = []
    for fget in getters:
        if fget in used_getters:
            continue
        for fset in setters:
            if fset in used_setters:
                continue
            property_ = self.recognizer.create_property( fget, fset )
            if property_:
                print "A PROPERTY!"
                if self.__is_legal_property( property_ ):
                    used_getters.add( fget )
                    used_setters.add( fset )
                    properties.append( property_ )
                    break
                else:
                    self.__report_illegal_property( property_ )
    return properties

  def __is_legal_property( self, property_ ):
    """property is legal if it does not hide other declarations"""
    return True

    print "CHECK LEGALITY!"

    def is_relevant( decl ):
      irrelevant_classes = ( declarations.constructor_t
                             , declarations.destructor_t
                             , declarations.typedef_t )
      if isinstance( decl, irrelevant_classes ):
          print "irrelevant class"
          return False
      if decl.ignore:
          print "decl.ignore"
          return False
      if decl.alias != property_.name:
          print "decl.alias != property_.name:"
          return False
      if self.exclude_accessors and ( decl is property_.fget or decl is property_.fset ):
          print "decl is property_.fget or decl is property_.fset ):"
          return False
      return True

    relevant_decls = []
    relevant_classes = [self.cls] + self.recognizer.base_classes( self.cls )
    for cls in relevant_classes:
        relevant_decls.extend( cls.decls( is_relevant, recursive=False, allow_empty=True ) )
    return not bool( relevant_decls )



def member_function_code_creators():
  return ( pyplusplus.code_creators.calldef.mem_fun_overloads_class_t,
           pyplusplus.code_creators.calldef.mem_fun_overloads_t,
           pyplusplus.code_creators.calldef.mem_fun_private_pv_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_private_v_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_pv_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_pv_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_s_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_s_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_v_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_v_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_protected_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_pv_t,
           pyplusplus.code_creators.calldef.mem_fun_pv_wrapper_t,
           pyplusplus.code_creators.calldef.mem_fun_t,
           pyplusplus.code_creators.calldef.mem_fun_v_t,
           pyplusplus.code_creators.calldef.mem_fun_v_wrapper_t
          )








