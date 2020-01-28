using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace Reconstruct
{
	struct ScaffIndices
	{
		public int st;
		public int ed;
		public bool reverse;
	};
	class ScaffParams 
	{
		public ScaffParams(bool hasSt, bool hasEd, Allele st, Allele ed)
		{
			startJoin = hasSt;
			endJoin = hasEd;
			stAllele = st;
			edAllele = ed;
		}
		public bool startJoin;
		public bool endJoin;
		public Allele stAllele;
		public Allele edAllele;
	};
	public class MetaScaffold
	{
		Scaffold _starter;
		List<Scaffold> _pathScaffs;
		List<ScaffIndices> _path;
		public bool _modeReverse;
		int _index;

		public MetaScaffold (Scaffold start, int index)
		{
			_index = index;
			_starter = start;
			_starter._inMeta = true;
			_path = new List<ScaffIndices> ();
		}

		public int GetLength()
		{
			int accu = 0;
			foreach (ScaffIndices si in _path) {
				accu += si.ed - si.st;
			}
			return accu;
		}

		public void WriteScaffolds(StreamWriter sw)
		{
			string title = ">MetaScaff" + _index.ToString () + "_" + GetLength ().ToString ();
			sw.WriteLine (title);
			StringBuilder stb = new StringBuilder ();

			if (_modeReverse) {
				int len = _pathScaffs.Count;
				for (int i = len - 1; i >= 0; i--) {
					int slen = _path [i].ed - _path [i].st;
					int st = _path [i].st;
					string toAdd = _pathScaffs[i]._sequence.Substring(st, slen);
					if (_path [i].reverse) {
						toAdd = RevComp (ref toAdd);
					}
					stb.Append (toAdd);
				}
			} else {
				int len = _pathScaffs.Count;
				for (int i = 0; i < len; i++) {
					int slen = _path [i].ed - _path [i].st;
					int st = _path [i].st;
					string toAdd = _pathScaffs[i]._sequence.Substring(st, slen);
					if (_path [i].reverse) {
						toAdd = RevComp (ref toAdd);
					}
					stb.Append (toAdd);
				}
			}
			string seq = stb.ToString ();
			Scaffold.WriteSequence (sw, ref seq);
		}

		public string RevComp(ref string toAdd)
		{
			char[] array = toAdd.ToCharArray ();
			Array.Reverse (array);

			int len = array.Length;
			for (int i = 0; i < len; i++) {
				if (array [i] == 'A') {
					array [i] = 'T';
				} else if (array [i] == 'T') {
					array [i] = 'A';
				} else if (array [i] == 'C') {
					array [i] = 'G';
				} else if (array [i] == 'G') {
					array [i] = 'C';
				} else if (array [i] == 'a') {
					array [i] = 't';
				} else if (array [i] == 't') {
					array [i] = 'a';
				} else if (array [i] == 'g') {
					array [i] = 'c';
				} else if (array [i] == 'c') {
					array [i] = 'g';
				}
			}

			string next = new string (array);

			return next;
		}

		public int GenerateOverlapPath()
		{
			List<Scaffold> firstPath = new List<Scaffold> ();
			List<ScaffParams> pathParams = new List<ScaffParams> ();
			_path = new List<ScaffIndices> ();
			_pathScaffs = new List<Scaffold> ();
			Dictionary<string, int> siTest = new Dictionary<string, int> ();
			Dictionary<string, int> pathTest = new Dictionary<string, int> ();

			int cnt = 0;

			Scaffold current = _starter;
			Allele connect;
			bool fromStart = current.OnlyStartAvail ();
			bool connectEnd = false;

			connect = current.GetNextConnection (true, connectEnd);
			pathTest.Add (_starter._name, 1);
			if (fromStart) {
				_modeReverse = true;
				pathParams.Add (new ScaffParams (true, false, connect, null));
			} else {
				_modeReverse = false;
				pathParams.Add (new ScaffParams (false, true, null, connect));
			}
			firstPath.Add (current);

			connectEnd = !ConnectsToStart (connect);

			bool leaveLoop = false;
			if (connect._target._inMeta) {
				leaveLoop = true;
			}

			while (!leaveLoop) {

				fromStart = !connectEnd;

				connectEnd = !ConnectsToStart (connect);
				firstPath.Add (connect._target);
				current = connect._target;
				current._inMeta = true;

				connect = current.GetNextConnection (false, connectEnd);
				if (connect == null || pathTest.ContainsKey(connect._target._name) || connect._target._inMeta) {
					leaveLoop = true;
					pathParams.Add (new ScaffParams(false, false, null, null));
				} else {
					pathTest.Add (connect._target._name, 1);
					if (fromStart) {
						pathParams.Add (new ScaffParams (true, false, connect, null));
					} else {
						pathParams.Add (new ScaffParams (false, true, null, connect));
					}
				}
			}

			int len = firstPath.Count;

			bool previousExists = false;
			ScaffParams previous = null;
			bool doReverse = false;
			bool PastReverse = false;

			for (int i = 0; i < len; i++) {
				
				ScaffParams scf = pathParams [i];
				Scaffold sc = firstPath [i];

				int stTrim = 0; int edTrim = sc._length;
				int pad = 0;

				if (i == len - 1 || !siTest.ContainsKey (firstPath [i + 1]._name)) {
					if (scf.startJoin) {
						connect = scf.stAllele;
						pad = ((connect._ed - connect._st) / 2);
						stTrim = connect._st + pad;
						edTrim = sc._length;
					} else if (scf.endJoin) {
						connect = scf.edAllele;
						pad = ((connect._ed - connect._st) / 2);
						stTrim = 0;
						edTrim = connect._st + ((connect._ed - connect._st) / 2);
					}
				}

				if (previousExists) {
					bool toSt = false;
					Allele prev = null;
					if (previous.endJoin) {
						prev = previous.edAllele;
					} else if (previous.startJoin) {
						prev = previous.stAllele;
					}

					toSt = ConnectsToStart (prev);
					int smallest = Math.Min (prev._tSt, prev._tEd);
					int largest = Math.Max (prev._tEd, prev._tSt);

					if (toSt) {
						stTrim = smallest + ((largest - smallest) / 2);
					} else {
						edTrim = smallest + ((largest - smallest) / 2);
					}
				}

				if ((PastReverse ^ connect._reversed) && i != 0) {
					doReverse = true;
					PastReverse = true;
				} else {
					doReverse = false;
					PastReverse = false;
				}

				ScaffIndices si = new ScaffIndices ();

				si.st = stTrim;
				si.ed = edTrim;
				si.reverse = doReverse;

				if (!siTest.ContainsKey (sc._name)) {
					_path.Add (si);
					siTest.Add (sc._name, 0);
					_pathScaffs.Add (sc);
				}

				previous = scf;
				previousExists = true;
			}
			cnt = _pathScaffs.Count - 1;

			if (cnt == 0) {
				_starter._inMeta = false;
			}
			return cnt;
		}

		protected bool ConnectsToStart(Allele a)
		{
			Scaffold sref = a._target;

			int smallest = Math.Min (a._tEd, a._tEd);
			int largest = Math.Max (a._tEd, a._tSt);

			if (smallest < sref._length - largest) {
				return true;
			} else {
				return false;
			}
		}
	}
}

