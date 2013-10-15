package Fasterq;
use Inline C => Config => LIBS => '-lz';
use Inline C => <<'CODE';

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define SANGER 1
#define ILLUMINA 2
#define SOLEXA 3

SV* fasterq_init(char* class, char* file) {
    gzFile fp;
    kseq_t *seq;
    fp = gzopen(file, "r");
    seq = kseq_init(fp);    
    SV* obj_ref = newSViv(0);
    SV* obj = newSVrv(obj_ref, class);
    sv_setiv(obj, (IV)seq);
    SvREADONLY_on(obj);
    return obj_ref;
}

int read_seq(SV* obj) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    int l = kseq_read(seq);
    return l;
}

int get_min_qual(SV* obj, int encoding, int length) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    char* qual = seq->qual.s;
    int minqual;
    int i;
	if ( SANGER == encoding ) {
		minqual = 126 - 33;
		for ( i = 0; i < length; i++ ) {	
			char rescaled = qual[i] - 33;
			if ( rescaled < minqual ) {
				minqual = rescaled;
			}
		}
	}
	else if ( ILLUMINA == encoding ) {	
		minqual = 126 - 64;
		for ( i = 0; i < length; i++ ) {
			char rescaled = qual[i] - 64;
			if ( rescaled < minqual ) {
				minqual = rescaled;
			}
		}	
	}
	else if ( SOLEXA == encoding ) {
		minqual = 126 - 59;
		for ( i = 0; i < length; i++ ) {
			char rescaled = qual[i] - 64;
			if ( rescaled < minqual ) {
				minqual = rescaled;
			}
		}	
	}
	return minqual;
}

int read_seq_min_qual(SV* obj, int minqual, double p, int encoding, HV* seen) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    int length;
	unsigned int iseed = (unsigned int)time(NULL);
  	srand (iseed);    
	while( ( length = kseq_read(seq) ) > 0 ) {
		if ( ! hv_exists(seen, seq->seq.s, 100) ) {
			double roll = ( 1.0 * rand() ) / ( 1.0 * RAND_MAX ); // roll the dice
			int qual = get_min_qual(obj,encoding,length);
			if ( qual >= minqual && roll <= p  ) {
				hv_store(seen, seq->seq.s, 100, NULL, 0);
				return length;
			}
		}
	}
	return 0;
}

double get_mean_qual(SV* obj, int encoding, int length) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    char* qual = seq->qual.s;
    double result = 0.0;
    int i;
	if ( SANGER == encoding ) {
		for ( i = 0; i < length; i++ ) {	
			result += ( ( 1.0 * qual[i] ) - 33.0 );
		}
	}
	else if ( ILLUMINA == encoding ) {
		for ( i = 0; i < length; i++ ) {
			result += ( ( 1.0 * qual[i] ) - 64.0 );
		}	
	}
	else if ( SOLEXA == encoding ) {
		for ( i = 0; i < length; i++ ) {
			result += ( ( 1.0 * qual[i] ) - 64.0 );
		}	
	}
	return result / ( 1.0 * length );
}

int read_seq_mean_qual(SV* obj, double meanqual, double p, int encoding, HV* seen) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    int length;
	unsigned int iseed = (unsigned int)time(NULL);
  	srand (iseed);    
	while( ( length = kseq_read(seq) ) > 0 ) {
		if ( ! hv_exists(seen, seq->seq.s, 100) ) {
			double roll = ( 1.0 * rand() ) / ( 1.0 * RAND_MAX ); // roll the dice
			double qual = get_mean_qual( obj, encoding, length );
			if ( qual >= meanqual ) {
				if ( roll <= p ) {
					// remember that we have returned this
					hv_store(seen, seq->seq.s, 100, NULL, 0);
					return length;
				}
			}			
		}
	}
	return 0;
}

int read_seq_with_name(SV* obj, char* name) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    int length;
	while( ( length = kseq_read(seq) ) > 0 ) {
		if ( strcmp(seq->name.s,name) == 0 ) {
			return length;
		}
	}
	return 0;
}

void finalize(SV* obj) {
    kseq_destroy((kseq_t*)SvIV(SvRV(obj)));
}

char* get_name(SV* obj) {
    return ((kseq_t*)SvIV(SvRV(obj)))->name.s;
}

char* get_seq(SV* obj) {
    return ((kseq_t*)SvIV(SvRV(obj)))->seq.s;   
}

char* get_qual(SV* obj) {
    return ((kseq_t*)SvIV(SvRV(obj)))->qual.s;    
}

char* get_comment(SV* obj) {
    return ((kseq_t*)SvIV(SvRV(obj)))->comment.s;        
}

int to_string(SV* obj,PerlIO* handle) {
    kseq_t* seq = (kseq_t*)SvIV(SvRV(obj));
    char result[50 + 100 + 1 + 100 + 4];
    snprintf(result, sizeof result, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);
    return PerlIO_puts(handle,result);
}

CODE

sub SANGER () { 1 }
sub ILLUMINA () { 2 }
sub SOLEXA () { 3 }

sub new {
    my $class = shift;
    my $file = shift;
    return $class->fasterq_init($file);
}

sub next_seq {
	my ( $self, %args ) = @_;
	my $encoding    = $args{'-encoding'} || ILLUMINA;
	my $probability = $args{'-probability'} || 1.00;
	if ( $args{'-meanqual'} ) {
		my $l = $self->read_seq_mean_qual($args{'-meanqual'},$probability,$encoding,$args{'-seen'}||{});
		return $l ? $self : undef;
	}
	elsif ( $args{'-minqual'} ) {
		my $l = $self->read_seq_min_qual($args{'-minqual'},$probability,$encoding,$args{'-seen'}||{});
		return $l ? $self : undef;
	}
	elsif ( $args{'-minlength'} ) {
		my $l = $self->read_seq_min_length($args{'-minlength'},$probability);
		return $l ? $self : undef;
	}
	elsif ( $args{'-name'} ) {
		my $l = $self->read_seq_with_name($args{'-name'});
		return $l ? $self : undef;
	}
}

sub DESTROY { shift->finalize }

1;
