#
#
#
class TableFormat:

  def __init__(self, dc, args={}):
      if not args.has_key('showmode'):
         if not args.has_key('columns'): raise dc.MissingArgument('table', 'columns')
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
          config[col] = seld.dc.condreducer(lspl[k])
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
      for col in self.columns: f.write('%s  ' % self.dc.condformatter(col))
      f.write('\n')

  def WriteConfig(self, f, config):
      for col in self.columns: 
          if not config.has_key(col): raise self.dc.MissingTag('table', col)
          f.write('%s  ' % (self.dc.condformatter(config[col])))
      f.write('\n')

  def WriteMany(self, f, configs):
      for c in configs: self.WriteConfig(f, c)

#
#
#
def pluginLoader(dc, args={}): return TableFormat(dc, args)

