//
// —”¶¬ƒNƒ‰ƒX by g0307 2006
//

#ifndef _RNG_
#define _RNG_	1

class RNG {

public:

	RNG(){}

	// [0,1]‚Ì—”‚ğ¶¬
	float random() const
	{
		return xor128()*(1.0/4294967296.0);
	}
	
	inline unsigned long xor128() const{
		static unsigned long x=123456789,y=362436069,z=521288629,w=88675123;
		unsigned long t;
		t=(x^(x<<11));x=y;y=z;z=w; return( w=(w^(w>>19))^(t^(t>>8)) );
	}
};

#endif
