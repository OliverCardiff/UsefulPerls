using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace Generate_Candidates
{
	public class Genome
	{
		public Dictionary<string, Scaffold> _scaffolds;
		string _gfile;
		string _wFile;
		string _bFile;

		public Genome (string file, string wFile, string bFile)
		{
			_gfile = file;
			_wFile = wFile;
			_bFile = bFile;
			_scaffolds = new Dictionary<string, Scaffold> ();
		}

		public void GenerateAllAlleles()
		{
			int cnt = 0;
			int increment = 50000000;
			int threshold = increment;
			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				if (cnt > threshold) {
					threshold += increment;
					Console.WriteLine ("Processed " + cnt + " bases in Allele generation");
				}

				kvp.Value.GenerateAlleles ();
				cnt += kvp.Value._length;
			}
		}

		public void WriteLinksFile(string name)
		{
			using (StreamWriter sw = new StreamWriter (name)) {
				foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
					kvp.Value.WriteLinks (sw);
				}
			}
		}

		public void GenerateAllAnchors()
		{
			int cnt = 0;
			int increment = 10000000;
			int threshold = increment;
			int anCnt = 0;

			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				if (cnt > threshold) {
					threshold += increment;
					Console.WriteLine ("Processed " + cnt + " bases, " + anCnt + " total spans");
				}
				anCnt += kvp.Value.PileupBlasts();
				cnt += kvp.Value._length;
			}
		}

		public bool LoadWigs()
		{
			string line = "";
			bool retVal = false;
			string[] spsSp;
			string[] spsEq;
			string[] spsTab;
			int cnt = 0;

			Scaffold current = null;
			List<int> curWig = new List<int> ();

			char[] spChar = new char[1] { ' ' };
			using (StreamReader sr = new StreamReader (_wFile)) {
				while ((line = sr.ReadLine ()) != null) {
					if (cnt > 0) {
						retVal = true;
						spsTab = line.Split ('\t');

						if (cnt % 10000000 == 0) {
							Console.WriteLine ("Read " + cnt + " lines of wig!");
						}

						if (spsTab.Length < 2 && line.Length > 5) {
								
							spsSp = line.Split (spChar, StringSplitOptions.RemoveEmptyEntries);
							spsEq = spsSp[1].Split ('=');
							if (current != null) {
								current.GiveWig (curWig);
							}
							if (_scaffolds.ContainsKey (spsEq [1])) {
								current = _scaffolds [spsEq [1]];
							} else {
								current = null;
							}

							curWig = new List<int> ();
						} else if (spsTab.Length == 2) {
							int cov;
							if (int.TryParse (spsTab [1], out cov)) {
								curWig.Add (cov);
							}
						}
					}
					cnt++;
				}
				if (current != null) {
					current.GiveWig (curWig);
				}
			}

			return retVal;
		}

		public bool LoadBlast()
		{
			string line = "";
			bool retVal = false;
			int cnt = 0;

			using (StreamReader sr = new StreamReader (_bFile)) {
				while ((line = sr.ReadLine ()) != null) {
					retVal = true;
					string[] spstr = line.Split ('\t');
					cnt++;

					if (cnt % 500000 == 0) {
						Console.WriteLine ("Read " + cnt + " lines of blast!");
					}

					bool parseSuccess = true;
					double pid;
					int st; int ed;
					int tst; int ted;

					//Starts and ends
					if (!int.TryParse (spstr [6], out st)) {
						parseSuccess = false;
					}
					if (!int.TryParse (spstr [7], out ed)) {
						parseSuccess = false;
					}
					if (!int.TryParse (spstr [8], out tst)) {
						parseSuccess = false;
					}
					if (!int.TryParse (spstr [9], out ted)) {
						parseSuccess = false;
					}

					if (!double.TryParse (spstr [2], out pid)) {
						parseSuccess = false;
					}

					if (parseSuccess && pid > MainClass.BLASTPERCID) {
						Scaffold Q = _scaffolds [spstr [0]];
						Scaffold T = _scaffolds [spstr [1]];

						Blast next = new Blast (Q, T, pid, st, ed, tst, ted);
						Q.GiveBlast (next);
					}
				}
			}
			return retVal;
		}

		public bool LoadGenome()
		{
			string line = "";
			string sequence = "";
			string id = "";
			bool retVal = false;
			int count = 0;
			int thresh = 50000000;
			int interval = 50000000;
			StringBuilder stb = new StringBuilder ();

			using (StreamReader sr = new StreamReader (_gfile)) {
				while ((line = sr.ReadLine ()) != null) {
					retVal = true;

					if (!line.StartsWith (">")) {
						stb.Append (line);
					} else {
						if (stb.Length > 0) {
							sequence = stb.ToString ();
							_scaffolds.Add (id, new Scaffold (id, ref sequence));
							count += sequence.Length;

							if (count > thresh) {
								thresh += interval;
								Console.WriteLine ("Bases read = " + count);
							}

							stb = new StringBuilder ();
						}
						id = line.TrimStart ('>');
					}
				}
				if (!_scaffolds.ContainsKey (id)) {
					_scaffolds.Add (id, new Scaffold (id, ref sequence));
				}
			}
			return retVal;
		}
	}
}

