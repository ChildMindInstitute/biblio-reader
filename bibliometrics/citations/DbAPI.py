import MySQLdb, MySQLdb.cursors 
class DbAPI:
	def __init__(self):
		self.con = MySQLdb.connect(host='127.0.0.1', user='root', passwd='renew', db='annotate', cursorclass=MySQLdb.cursors.DictCursor, connect_timeout=3, charset='utf8')
		self.con.autocommit(True)
		self.cur = self.con.cursor()		

	def execute(self,sqlString):
		if (sqlString!=''):
			self.cur.execute(sqlString)
			return self.cur.fetchall()
		else:
			return False

	def singleRowQuery(self,sqlString):
		if (sqlString!=''):
			self.cur.execute(sqlString)
			return self.cur.fetchone()
		else:
			return False
		
	def close(self):
		return self.con.close()

	def commit(self):
		return self.con.commit()
