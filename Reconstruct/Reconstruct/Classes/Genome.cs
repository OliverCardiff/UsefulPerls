using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace Reconstruct
{
	public class Genome
	{
		public Dictionary<string, Scaffold> _scaffolds;
		public List<Scaffold> _singleJoins;
		public List<MetaScaffold> _metaScaffs;
		string _gfile;
		string _lFile;
		string _outFile;

		public Genome (string file, string linksFile, string outFile)
		{
			_gfile = file;
			_lFile = linksFile;
			_outFile = outFile;
			_scaffolds = new Dictionary<string, Scaffold> ();
			_singleJoins = new List<Scaffold> ();
			_metaScaffs = new List<MetaScaffold> ();
		}

		public void WriteAlleles()
		{
			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				kvp.Value.WriteAllAlleles ();
			}	
		}

		public void BuildMetaScaffolds()
		{
			int cnt = 0;
			int joins = 0;
			int maxo = 0;
			int len = _singleJoins.Count;
			for (int i = 0; i < len; i++) {
				if (_singleJoins [i]._inMeta == false) {
					cnt++;
					MetaScaffold meta = new MetaScaffold (_singleJoins [i], cnt);
					int join = meta.GenerateOverlapPath ();
					joins += join;
					if (maxo < join + 1) {
						maxo = join + 1;
					}
					if (join > 0) {
						_metaScaffs.Add (meta);
					}
					/*if (cnt % 10 == 0) {
						Console.WriteLine ("Made " + cnt + " metas so far");
					}*/
				}
			}
			Console.WriteLine ("Made MetaScaffolds: " + cnt);
			Console.WriteLine ("Made Total Joins: " + joins);
			Console.WriteLine ("Highest Meta Scaff Count: " + maxo + "\n");
		}

		public void KeepJoins()
		{
			int cnt = 0;
			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				int joins = kvp.Value.KeepOnlyJoins ();
				cnt += joins;
				if (joins == 1) {
					_singleJoins.Add (kvp.Value);
				}
			}
			_singleJoins.Sort ();

			Console.WriteLine ("All Confirmed Join Count: " +  cnt.ToString() + "\n");
			Console.WriteLine ("Single Join Scaffolds: " +  _singleJoins.Count.ToString() + "\n");
		}

		public void ExpandReflections(bool removeContained)
		{
			int cnt = 0;
			int alCnt = 0;
			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				alCnt += kvp.Value.AlleleCount;
				cnt += kvp.Value.ReflectionExpand ();
			}

			Console.WriteLine ("Non-reflected links removed: " + (alCnt - cnt).ToString());
			Console.WriteLine ("Number of links 'Reflection-expanded': " + cnt.ToString());

			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				kvp.Value.SwitchAllelesToKeeps ();
			}

			if (removeContained) {
				int cnt2 = 0;
				foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
					cnt2 += kvp.Value.RemoveReflectionContained ();
				}

				Console.WriteLine ("Reflection mismatches removed: " + cnt2 + "\n");

				foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
					kvp.Value.SwitchAllelesToKeeps ();
				}
			} else {
				Console.WriteLine ();
			}
		}

		public void MergeOverlaps()
		{
			int cnt = 0;
			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				cnt += kvp.Value.MergeOverlaps ();
			}

			Console.WriteLine ("Overlapping Spans Merged: " + cnt + "\n");
		}

		public void EliminateSmallFrags()
		{
			List<string> toremove = new List<string> ();

			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				if (kvp.Value.AssessForRemoval ()) {
					toremove.Add (kvp.Key);
				}
			}
			foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
				kvp.Value.RemoveKilledLinks (toremove);
			}

			foreach (string st in toremove) {
				_scaffolds.Remove (st);
			}

			Console.WriteLine ("Removed " + toremove.Count + " small frags\n");
		}

		public bool LoadAlleles()
		{
			string line = "";
			bool retVal = false;

			using (StreamReader sr = new StreamReader (_lFile)) {
				while ((line = sr.ReadLine ()) != null) {
					retVal = true;
					string[] spstr = line.Split ('\t');

					bool parseSuccess = true;
					int st; int ed;
					int tst; int ted;

					//Starts and ends
					if (!int.TryParse (spstr [2], out st)) {
						parseSuccess = false;
					}
					if (!int.TryParse (spstr [3], out ed)) {
						parseSuccess = false;
					}
					if (!int.TryParse (spstr [4], out tst)) {
						parseSuccess = false;
					}
					if (!int.TryParse (spstr [5], out ted)) {
						parseSuccess = false;
					}

					if (parseSuccess) {
						int l1 = ed - st;
						int l2 = Math.Abs (tst - ted);
						if (l1 > 100 && l2 > 100) {
							Scaffold Q = _scaffolds [spstr [0]];
							Scaffold T = _scaffolds [spstr [1]];

							Allele next = new Allele (Q, T, st, ed, tst, ted);
							Q.GiveAllele (next);
						}
					}
				}
			}
			return retVal;
		}

		public bool WriteGenome(bool actuallyDoIt)
		{
			bool retVal = false;
			int scCnt = 0;
			int mtCnt = 0;
			int outLen = 0;
			int inLen = 0;

			using (StreamWriter sw = new StreamWriter (_outFile)) {
				if (sw != null) {
					retVal = true;
				}
				foreach (MetaScaffold ms in _metaScaffs) {
					if(actuallyDoIt) {
						ms.WriteScaffolds (sw);
					}
					inLen += ms.GetLength ();
					mtCnt++;
				}
				foreach (KeyValuePair<string, Scaffold> kvp in _scaffolds) {
					if (!kvp.Value._inMeta) {
						if(actuallyDoIt) {
							kvp.Value.WriteOwnSeq (sw);
						}
						outLen += kvp.Value._length;
						scCnt++;
					}
				}
			}
			Console.WriteLine ("\nMeta Count:\t" + mtCnt);
			Console.WriteLine ("Meta Bases:\t" + inLen);
			Console.WriteLine ("Scaff Count:\t" + scCnt);
			Console.WriteLine ("Scaff Bases:\t" + outLen);

			return retVal;
		}

		public bool LoadGenome()
		{
			string line = "";
			string sequence = "";
			string id = "";
			bool retVal = false;
			int count = 0;
			int interval = 100000000;
			int thresh = interval;
			StringBuilder stb = new StringBuilder ();

			using (StreamReader sr = new StreamReader (_gfile)) {
				while ((line = sr.ReadLine ()) != null) {
					retVal = true;

					if (!line.StartsWith (">")) {
						stb.Append (line);
					} else {
						if (stb.Length > 0) {
							sequence = stb.ToString ();
							if (_scaffolds.ContainsKey (id)) {
								_scaffolds[id].GiveSeq(ref sequence);
							}
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
				if (_scaffolds.ContainsKey (id)) {
					_scaffolds[id].GiveSeq(ref sequence);
				}
			}
			return retVal;
		}

		public bool LiteLoadGenome()
		{
			string line = "";
			string id = "";
			bool retVal = false;

			using (StreamReader sr = new StreamReader (_gfile + ".fai")) {
				while ((line = sr.ReadLine ()) != null) {
					retVal = true;
					bool parseSuccess = true;
					int len = 0;
					string[] spstr = line.Split ('\t');

					if (!int.TryParse (spstr [1], out len)) {
						parseSuccess = false;
					}
					id = spstr [0];

					if (!_scaffolds.ContainsKey (id) && parseSuccess) {
						_scaffolds.Add (id, new Scaffold (id, len));
					}
				}
			}
			return retVal;
		}
	}
}

