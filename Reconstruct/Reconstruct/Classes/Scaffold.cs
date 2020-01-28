using System;
using System.Collections.Generic;
using System.IO;

namespace Reconstruct
{
	/*public struct MergeUpdate
	{
		public Scaffold A;
		public Scaffold B;
		public Scaffold Next;
		Allele JoinBy;

		public MergeUpdate(Scaffold a, Scaffold b, Scaffold final, Allele joiner)
		{
			A = a;
			B = b;
			Next = final;
			JoinBy = joiner;
		}

	};*/
	public class Scaffold : IComparable
	{
		public override string ToString ()
		{
			return _name;
		}
		protected List<Allele> _alleles;
		protected List<Allele> _tokeep;

		public List<Allele> _startJoins;
		public List<Allele> _endJoins;

		public string _name;
		public int _length;
		public string _sequence;
		public bool _toKill;
		public bool _inMeta;
		public bool _mandatory;

		public int AlleleCount { get { return _alleles.Count; } }

		#region IComparable implementation

		public int CompareTo (object obj)
		{
			Scaffold oth = (Scaffold)obj;

			if (oth._length < this._length) {
				return -1;
			} else if (oth._length > this._length) {
				return 1;
			} else {
				return 0;
			}

			return 0; 
		}

		#endregion

		public Scaffold(string id, ref string seq)
		{
			_sequence = seq;
			_name = id;
			_length = seq.Length;
			_alleles = new List<Allele> ();
			_tokeep = new List<Allele> ();
			_inMeta = false;
			_mandatory = false;
		}

		public Scaffold(string id, int length)
		{
			_sequence = "";
			_name = id;
			_length = length;
			_alleles = new List<Allele> ();
			_tokeep = new List<Allele> ();
			_inMeta = false;
			_mandatory = false;
		}

		public bool OnlyStartAvail()
		{
			if (_startJoins.Count > 0) {
				return true;
			} else {
				return false;
			}
		}

		public Allele GetNextConnection(bool AnySingle, bool StartAvailable)
		{
			if (AnySingle) {
				if (_startJoins.Count > 0) {
					return _startJoins [0];
				} else if (_endJoins.Count > 0) {
					return _endJoins [0];
				} else {
					return null;
				}
			} else if (StartAvailable) {
				if (_startJoins.Count > 0) {
					_startJoins.Sort ();
					return _startJoins [0];
				} else {
					return null;
				}
			} else {
				if (_endJoins.Count > 0) {
					_endJoins.Sort ();
					return _endJoins [0];
				} else {
					return null;
				}
			}
		}

		public void WriteOwnSeq(StreamWriter sw)
		{
			sw.WriteLine (">" + _name);
			WriteSequence (sw, ref _sequence);
		}

		public static void WriteSequence(StreamWriter sw, ref string sequence)
		{
			int len = sequence.Length;
			int head = 0;
			while (head < len) {
				if (head + 100 < len) {
					sw.WriteLine (sequence.Substring (head, 100));
				} else {
					int linLen = len - head;
					sw.WriteLine (sequence.Substring (head, linLen));
				}
				head += 100;
			}
		}

		public void GiveSeq(ref string sequence)
		{
			_sequence = sequence;
		}

		public void GiveAllele(Allele a)
		{
			_alleles.Add (a);
		}

		public void SwitchAllelesToKeeps()
		{
			_alleles = _tokeep;
			_tokeep = new List<Allele> ();
		}

		public int KeepOnlyJoins()
		{
			int retVal = 0;
			_startJoins = new List<Allele> ();
			_endJoins = new List<Allele> ();

			foreach (Allele al in _alleles) {
				Scaffold tar = al._target;
				int stEdge = al._st;
				int edEdge = _length - al._ed;
				int tstEdge = al._tSt;
				int tedEdge = tar._length - al._tEd;

				if (al._reversed) {
					tstEdge = tar._length - al._tSt;
					tedEdge = al._tEd;
				}

				if (al._st < _length / MainClass.LengthDiv) {

					if (stEdge < tstEdge && tedEdge < edEdge) {
						_startJoins.Add (al);
						retVal++;
					}

				} else if  (al._ed > _length - (_length / MainClass.LengthDiv)) {
					if (edEdge < tedEdge && tstEdge < stEdge) {
						_endJoins.Add (al);
						retVal++;
					}
				}
			}

			return retVal;
		}

		public int RemoveReflectionContained()
		{
			int retVal = 0;

			foreach (Allele al in _alleles) {
				Scaffold t = al._target;
				if (!t.ReflectionMismatch(al)) {
					_tokeep.Add (al);
				} else {
					retVal++;
				}
			}

			return retVal;
		}

		public int MergeOverlaps()
		{
			int retVal = 0;

			int len = _alleles.Count;
			List<Allele> final = new List<Allele> ();

			for (int i = 0; i < len - 1; i++) {
				for (int j = i + 1; j < len; j++) {
					int prox = _alleles [i].GetProximity (_alleles [j]);
					int lenBoth = Math.Max(_alleles[i].Length, _alleles[j].Length);
					int bound = Math.Min (lenBoth / 10, 2000);

					if (_alleles[i]._target._name.Equals(_alleles [j]._target._name) && (_alleles [i].Overlaps(_alleles[j]) || prox < bound)) {
						_alleles [j] = new Allele (_alleles [i], _alleles [j]);
						retVal++;
						_alleles [i] = null;
						break;
					}
				}
			}

			foreach (Allele a in _alleles) {
				if (a != null) {
					final.Add (a);
				}
			}

			_alleles = final;

			return retVal;
		}

		public int ReflectionExpand()
		{
			int retVal = 0;

			foreach (Allele al in _alleles) {
				Scaffold t = al._target;
				List<Allele> matches = t.HasCloseMatch (al);
				if (matches.Count > 0) {
					retVal++;
					int smallEnd = al._st;
					int bigEnd = al._ed;
					int smallRef = al._tSt;
					int bigRef = al._tEd;

					foreach (Allele m in matches) {
						int tst = m._tSt; int refSt = m._st;
						int ted = m._tEd; int refEd = m._ed;
						if (m._reversed) {
							tst = m._tEd;
							ted = m._tSt;
						}
						if (al._reversed) {
							smallRef = Math.Max (refEd, smallRef);
							bigRef = Math.Min (bigRef, refSt);
						} else {
							smallRef = Math.Min (refSt, smallRef);
							bigRef = Math.Max (bigRef, refEd);
						}

						smallEnd = Math.Min (smallEnd, tst);
						bigEnd = Math.Max (bigEnd, ted);
					}
					al._st = smallEnd;
					al._ed = bigEnd;
					al._tSt = smallRef;
					al._tEd = bigRef;
					_tokeep.Add (al);
				}
			}
			return retVal;
		}

		protected bool ReflectionMismatch(Allele a)
		{
			int testSt = 0; int testEd = 0;
			string tnam = a._ref._name;
			if (a._reversed) {
				testSt = a._tSt;
				testEd = a._tEd;
			} else {
				testSt = a._tSt;
				testEd = a._tEd;
			}

			foreach (Allele al in _alleles) {
				if (al.Contains (testSt, testEd)) {
					if (!al._target.Equals (tnam)) {
						return true;
					}
				}
			}
			return false;
		}

		public void WriteAllAlleles()
		{
			foreach (Allele a in _alleles) {
				a.WriteToScreen ();
			}
		}

		public List<Allele> HasCloseMatch(Allele a)
		{
			int testSt = 0; int testEd = 0;
			string tnam = a._ref._name;
			if (a._reversed) {
				testSt = a._tSt;
				testEd = a._tEd;
			} else {
				testSt = a._tSt;
				testEd = a._tEd;
			}

			double boundary = Math.Max(3000, Math.Min (((double)a.Length), 10000));

			List<Allele> retVal = new List<Allele> ();

			foreach (Allele al in _alleles) {
				if (al.Overlaps (testSt, testEd)) {
					if (al._target._name.Equals (tnam)) {
						retVal.Add (al);
					}
				} else if ((double)al.GetProximity(testSt, testEd) < boundary) {
					if (al._target._name.Equals (tnam)) {
						retVal.Add (al);
					}
				} 
			}
			return retVal;
		}

		public bool AssessForRemoval()
		{
			double test_len = (double)_length * 0.8;
			double thresh_len = (double)_length * 0.35;

			if (_alleles.Count == 1) {
				double len2 = (double)(_alleles [0].Length);

				if (len2 > test_len && _alleles[0]._target._length > _length) {
					return true;
				} else if (len2 > thresh_len && IsSurrounded (_alleles [0])) {
					return true;
				}
			} 
			return false;
		}

		protected bool IsSurrounded(Allele al)
		{
			int stdiff = al._st;
			int eddiff = _length - al._ed;

			int targetSt = 0; int targetEd = 0;

			if (!al._reversed) {
				targetSt = al._tSt;
				targetEd = al._target._length - al._tEd;
			} else {
				targetSt = al._target._length - al._tSt;
				targetEd = al._tEd;
			}

			if (stdiff < targetSt && eddiff < targetEd) {
				return true;
			} else {
				return false;
			}
		}

		public void RemoveKilledLinks(List<string> toKill)
		{
			List<Allele> survivors = new List<Allele> ();
			foreach (Allele al in _alleles) {
				bool keep = true;

				foreach (string st in toKill) {
					if (al._target._name.Equals (st)) {
						keep = false;
						break;
					}
				}
				if (keep) {
					survivors.Add (al);
				}
			}

			_alleles = survivors;
		}

		protected List<Allele> JoinAdjacentAlleles(List<Allele> all)
		{
			int len = all.Count;
			List<Allele> final = new List<Allele> ();

			for (int i = 0; i < len - 1; i++) {
				for (int j = i + 1; j < len; j++) {
					int prox = all [i].GetProximity (all [j]);
					int lenBoth = Math.Max(all[i].Length, all[j].Length);
					int bound = Math.Min (lenBoth / 10, 2000);

					if (prox < bound) {
						all [j] = new Allele (all [i], all [j]);
						all [i] = null;
						break;
					}
				}
			}

			foreach (Allele a in all) {
				if (a != null) {
					final.Add (a);
				}
			}

			return final;
		}

		protected List<Allele> MergeAlleleLists(List<Allele> broad, List<Allele> narrow)
		{
			int nLen = narrow.Count;
			int bLen = broad.Count;

			for (int j = 0; j < bLen; j++) 
			{
				for (int i = 0; i < nLen; i++) {
					if (narrow [i] != null) {
						if (broad [j].Contains (narrow [i])) {
							narrow [i] = null;
						} else if (broad [j].Overlaps (narrow [i])) {
							broad [j] = new Allele (broad [j], narrow [i]);
							narrow [i] = null;
						}
					}
				}
			}

			return broad;
		}
	}
}

