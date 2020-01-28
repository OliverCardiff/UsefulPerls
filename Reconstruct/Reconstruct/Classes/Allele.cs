using System;

namespace Reconstruct
{
	public class Allele : IComparable
	{
		public Scaffold _ref;
		public Scaffold _target;

		public int _st;
		public int _ed;
		public int _tSt;
		public int _tEd;
		public bool _reversed;

		public int Length { get { return _ed - _st; } }
		public int Oth_Length { get {
				if (_reversed) {
					return _tSt - _tEd;
				} else {
					return _tEd - _tSt;
				}
				} }
	
		public Allele (Scaffold sref, Scaffold oref, int st, int ed, int tst, int ted)
		{
			_ref = sref;
			_target = oref;
			_st = st;
			_ed = ed;
			_tEd = ted;
			_tSt = tst;
			if (_tSt > _tEd) {
				_reversed = true;
			} else {
				_reversed = false;
			}
		}

		#region IComparable implementation

		public int CompareTo (object obj)
		{
			Allele oth = (Allele)obj;

			if (oth.Length < this.Length) {
				return -1;
			} else if (oth.Length > this.Length) {
				return 1;
			} else {
				return 0;
			}

			return 0; 
		}

		#endregion

		public void WriteToScreen()
		{
			Console.WriteLine (_ref._name + "\t" + _target._name + "\t" + _st + "\t" + _ed + "\t" + _tSt + "\t" + _tEd);
		}

		public Allele(Allele A, Allele B)
		{
			if (A._reversed) {
				_tSt = Math.Max (A._tSt, B._tSt);
				_tEd = Math.Min (A._tEd, B._tEd);
			} else {
				_tSt = Math.Min (A._tSt, B._tSt);
				_tEd = Math.Max (A._tEd, B._tEd);
			}
			_st = Math.Min (B._st, A._st);
			_ed = Math.Max (B._ed, A._ed);
			_ref = A._ref;
			_target = A._target;
		}
			

		protected void ExpandToStart(Allele asp)
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

		protected void ExpandToEnd(Allele asp)
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

		public bool Overlaps(Allele oth)
		{
			if ((_st <= oth._ed && _st >= oth._st) || (_ed <= oth._ed && _ed >= oth._st)) {
				return true;
			} else {
				return false;
			}
		}

		public bool Overlaps(int othSt, int othEd)
		{
			if ((_st <= othEd && _st >= othSt) || (_ed <= othEd && _ed >= othSt)) {
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

		public bool Contains(int othSt, int othEd)
		{
			if (_st >= othSt && _ed <= othEd) {
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

		public int GetProximity(int othSt, int othEd)
		{
			int testA = othSt - _ed;
			int testB = _st - othEd;

			return Math.Max (testA, testB);
		}
	}
}

