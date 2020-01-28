using System;
using System.Collections.Generic;
using System.IO;

namespace Generate_Candidates
{
	public class Allele
	{
		protected List<AnchorBase> _bases;
		protected List<Blast> _blasts;
		protected Scaffold _ref;

		public int _st;
		public int _ed;

		public int Length { get { return _ed - _st; } }
		public List<AnchorSpan> _spans;

		public Allele (Scaffold sref, int st, int ed)
		{
			_ref = sref;
			_st = st;
			_ed = ed;
			_bases = new List<AnchorBase> ();
			_blasts = null;
		}

		public Allele(Allele Broad, Allele Narrow)
		{
			_st = Math.Min (Broad._st, Narrow._st);
			_ed = Math.Max (Broad._ed, Narrow._ed);
			_ref = Broad._ref;
			_bases = new List<AnchorBase> ();
			_blasts = null;
		}

		public void AddBlast(Blast b)
		{
			if (_blasts == null) {
				_blasts = new List<Blast> ();
			}
			_blasts.Add (b);
		}

		public void WriteSpans(StreamWriter sw)
		{
			foreach (AnchorSpan asp in _spans) {
				asp.Write (sw, _ref);
			}
		}

		public int GenerateAnchorBases()
		{
			int spanInt = 0;
			_spans = new List<AnchorSpan> ();

			for (int i = _st; i <= _ed; i++) {
				_bases.Add (new AnchorBase (null, i));
			}

			if (_blasts != null) {
				foreach (AnchorBase ab in _bases) {
					foreach (Blast b in _blasts) {
						if (ab.IsWithin (b)) {
							if (ab._occupied && ab._isAnchor) {
								if (!ab._target._name.Equals (b._B._name, StringComparison.Ordinal)) {
									ab._isAnchor = false;
									ab._target = null;
								}
							} else if(!ab._occupied){
								ab.SetTargetBase (b);
								ab._isAnchor = true;
								ab._occupied = true;
								ab._target = b._B;
							}
						}
					}
				}

				Dictionary<string, List<AnchorBase>> spanData = new Dictionary<string, List<AnchorBase>> ();

				foreach (AnchorBase ab in _bases) {
					if (ab._isAnchor) {
						if (spanData.ContainsKey (ab._target._name)) {
							spanData [ab._target._name].Add (ab);
						} else {
							spanData.Add (ab._target._name, new List<AnchorBase> ());
							spanData [ab._target._name].Add (ab);
						}
					}
				}

				List<AnchorSpan> spansA = new List<AnchorSpan> ();

				foreach (KeyValuePair<string, List<AnchorBase>> kvp in spanData) {
					AnchorSpan asp = new AnchorSpan (kvp.Value);
					if (asp.Length > 40 && asp.Count > 20) {
						spansA.Add (asp);
					}
				}

				if (spansA.Count > 1) {
					_spans = ProcessOverlaps (spansA);
				} else {
					_spans = spansA;
				}
				ExpandSpans (_spans);
				foreach (AnchorSpan aS in _spans) {
					spanInt++;
				}
			}
			return spanInt;
		}

		protected void ExpandSpans(List<AnchorSpan> spn)
		{
			if (spn.Count != 0) {
				if (spn.Count == 1) {
					ExpandToEnd (spn [0]);
					ExpandToStart (spn [0]);
				} else {

					int len = spn.Count;
					int tost = 0;
					int toed = 0;
					int min = 99999999;
					int max = 0;
					for (int i = 0; i < len; i++) {
						if (spn [i]._st < min) {
							min = spn [i]._st;
							tost = i;
						}
						if (spn [i]._ed > max) {
							max = spn [i]._ed;
							toed = i;
						}
					}

					ExpandToEnd (spn [toed]);
					ExpandToStart (spn [tost]);
				}
			}
		}

		protected void ExpandToStart(AnchorSpan asp)
		{
			if (asp._reversed == false) {
				int diff = asp._st - _st;
				if (diff < asp._tSt) {
					asp._st = _st;
					asp._tSt -= diff;
				} else {
					asp._st -= (asp._tSt - 1);
					asp._tSt = 1;
				}
			} else {
				int diff = asp._st - _st;
				int tdiff = asp._target._length - asp._tSt;

				if (diff < tdiff) {
					asp._st = _st;
					asp._tSt += diff;
				} else {
					asp._st -= tdiff;
					asp._tSt = asp._target._length;
				}
			}
		}

		protected void ExpandToEnd(AnchorSpan asp)
		{
			if (asp._reversed == false) {
				int diff = _ed - asp._ed;
				int tdiff = asp._target._length - asp._tEd;

				if (diff < tdiff) {
					asp._ed = _ed;
					asp._tEd += diff;
				} else {
					asp._ed += tdiff;
					asp._tEd = asp._target._length;
				}
			} else {
				int diff = _ed - asp._ed;

				if (diff < asp._tEd) {
					asp._ed = _ed;
					asp._tEd -= diff;
				} else {
					asp._ed += (asp._tEd - 1);
					asp._tEd = 1;
				}
			}
		}

		protected List<AnchorSpan> ProcessOverlaps(List<AnchorSpan> spn)
		{
			List<AnchorSpan> retVal = new List<AnchorSpan> ();
			List<AnchorSpan> sp2 = new List<AnchorSpan> ();
			List<bool> keeps = new List<bool> ();
			List<bool> kills = new List<bool> ();

			int len = spn.Count;

			for (int i = 0; i < len; i++) {
				keeps.Add (true);
			}

			for (int i = 0; i < len - 1; i++) {
				for (int j = i + 1; j < len; j++) {
					if (spn [j].Contains (spn [i])) {
						keeps [i] = false;
					} else if (spn [i].Contains (spn [j])) {
						keeps [j] = false;
					}
				}
			}

			for (int i = 0; i < len; i++) {
				if (keeps [i]) {
					sp2.Add (spn [i]);
					kills.Add (false);
				}
			}

			len = sp2.Count;

			if (len == 1) {
				retVal = sp2;
			} else {

				for (int i = 0; i < len - 1; i++) {
					for (int j = i + 1; j < len; j++) {
						if (kills [j] == false) {
							if (sp2 [i].Overlaps (sp2 [j])) {
								int ovSt;
								int ovEd;
								if (sp2 [i]._st < sp2 [j]._st) {
									ovSt = sp2 [j]._st;
									ovEd = sp2 [i]._ed;
								} else {
									ovSt = sp2 [i]._st;
									ovEd = sp2 [j]._ed;
								}
								
								int overlapLen = ovEd - ovSt;

								if (sp2 [i].Length - overlapLen < 300 ^ sp2 [j].Length - overlapLen < 300) {
									if (sp2 [i].Length < sp2 [j].Length) {
										kills [i] = true;
									} else {
										kills [j] = true;
									}
								} else if (sp2 [i].Length - overlapLen < 300 && sp2 [j].Length - overlapLen < 300) {
									double dI = sp2 [i].GetDensitySpan (sp2 [i]._st, sp2 [i]._ed);
									double dJ = sp2 [j].GetDensitySpan (sp2 [j]._st, sp2 [j]._ed);

									if (dI > dJ) {
										kills [j] = true;
									} else {
										kills [i] = true;
									}
								} else {
									double di = sp2 [i].GetDensitySpan (ovSt, ovEd);
									double dj = sp2 [j].GetDensitySpan (ovSt, ovEd);

									if (di > dj) {
										sp2 [j].DestroyOverlap (ovSt, ovEd);
									} else {
										sp2 [i].DestroyOverlap (ovSt, ovEd);
									}
								}
							}
						}
					}

				}

				for (int i = 0; i < len; i++) {
					if (!kills [i]) {
						retVal.Add (sp2 [i]);
					}
				}
			}

			return retVal;
		}

		public bool Overlaps(Allele oth)
		{
			if ((_st <= oth._ed && _st >= oth._st) || (_ed <= oth._ed && _ed >= oth._st)) {
				return true;
			} else {
				return false;
			}
		}

		public bool Contains(Allele oth)
		{
			if (_st >= oth._st && _ed <= oth._ed) {
				return true;
			} else {
				return false;
			}
		}

		public int GetProximity(Allele oth)
		{
			int testA = oth._st - _ed;
			int testB = _st - oth._ed;

			return Math.Max (testA, testB);
		}
	}
}

