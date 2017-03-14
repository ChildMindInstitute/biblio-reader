import os, sys, threading, subprocess

inDir = "/home/files/pdf/"
outDir = "/home/files/txt/"

class Command(object):
	def __init__(self, cmd):
		self.cmd = cmd
		self.process = None

	def run(self, timeout):
		def target():
			#print 'Thread started'
			self.process = subprocess.Popen(self.cmd, shell=True)
			self.process.communicate()
			#print 'Thread finished'

		thread = threading.Thread(target=target)
		thread.start()

		thread.join(timeout)
		if thread.is_alive():
			print 'Terminating process'
			self.process.terminate()
			thread.join()
		return self.process.returncode
		
def walkAndText(rootDir, outDir):
	for root, subFolders, files in os.walk(rootDir):
		for file in files:
			f = os.path.join(root,file)
			if f[-3:] == "pdf":
				print f
				outf = f
				outf = outf.replace(rootDir, "")
				outf = outf.replace("/", ".")
				outf = outf + ".txt"
				path = outDir + outf
				if os.path.exists(path) and os.path.getsize(path) > 0:
					print "Skipping", outf
				else:
					command = Command("pdftotext %s %s" % (f, path))
					command.run(10)
						  
if __name__ == "__main__":
	walkAndText(inDir, outDir)
 

