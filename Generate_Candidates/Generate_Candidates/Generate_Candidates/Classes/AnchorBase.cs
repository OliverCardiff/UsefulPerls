using System;

namespace Generate_Candidates
{
	public class AnchorBase
	{
		public bool _occupied;
		public bool _isAnchor;
		public Scaffold _target;
		public int _loc;
		public int _tarLoc;
		public bool _reversed;

		public AnchorBase (Scaffold target, int location)
		{
			_loc = location;
			_reversed = false;
			_target = target;
			_occupied = false;
			_isAnchor = false;
		}

		public bool IsWithin(Blast b)
		{
			if (_loc <= b._ed && _loc >= b._st) {
				return true;
			} else {
				return false;
			}
		}

		public bool IsWithin(int st, int ed)
		{
			if (_loc <= ed && _loc >= st) {
				return true;
			} else {
				return false;
			}
		}

		public void SetTargetBase(Blast b)
		{
			int inset = _loc - b._st;
			if (b._tEd > b._tSt) {
				_tarLoc = b._tSt + inset;
			} else {
				_tarLoc = b._tSt - inset;
				_reversed = true;
			}

		}
	}
}

