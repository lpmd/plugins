#
#
#

def condreducer(x):
    try: return int(x)
    except:
      try: return float(x)
      except: return x

def condformatter(x):
    try: return ('%8.8f' % x).ljust(12)
    except: 
      try: return ('%d' % x).ljust(12)
      except: return ('%s' % x).ljust(12)

class MissingArgument(Exception):

  def __init__(self, plugname, misarg):
      self.plugin = plugname
      self.arg = misarg
  
  def __str__(self): return ("[Error] (in plugin %s): Missing argument: \"%s\"" % (self.plugin, self.arg))

class TableFormat:

  def __init__(self, args={}):
      if not args.has_key('columns'): raise MissingArgument('table', 'columns')
      self.columns = args['columns'].split('/')

  def Show(self): print "*** table plugin information ***"

  def ReadHeader(self, f): pass          # there is no header to read

  def ReadConfig(self, f):
      config = {'NATOMS': 0}
      line = f.readline()
      if line == '': return None
      while line.strip().startswith('#'): line = f.readline()
      lspl, k = line.split(), 0
      for col in self.columns: 
          config[col] = condreducer(lspl[k])
          k = k + 1
      return config

  def ReadMany(self, f):
      data = []
      while True:
        c = self.ReadConfig(f)
        if c == None: break
        data.append(c)
      return data

  def WriteHeader(self, f):
      f.write('# ')
      for col in self.columns: f.write('%s  ' % condformatter(col))
      f.write('\n')

  def WriteConfig(self, f, config):
      for col in self.columns: f.write('%s  ' % (condformatter(config[col])))
      f.write('\n')

  def WriteMany(self, f, configs):
      for c in configs: self.WriteConfig(f, c)

#
#
#

def pluginLoader(args={}): return TableFormat(args)

