MODULE readsr_utils
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: input_string_len = 80

  PUBLIC :: keyword_contains
  PUBLIC :: index_of_delimiter
  PUBLIC :: readsr
  PUBLIC :: xstring
  PUBLIC :: readsi

CONTAINS

  ! ==================================================================
  FUNCTION keyword_contains(line,keyword,and,alias,but_not,cut_at) &
  RESULT (is_present)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), INTENT(in)         :: line, keyword
    CHARACTER(len=*), OPTIONAL, &
                      INTENT(in)         :: and, alias, but_not, cut_at

    CHARACTER(len=input_string_len)      :: next_words, copy_cut
 
    LOGICAL                              :: keyword_is_present, alias_is_present, &
                                            second_word_is_present, is_present, &
                                            exclusion_is_present
    INTEGER                              :: first, last

    keyword_is_present     = .false.
    second_word_is_present = .false.
    alias_is_present       = .false.
    exclusion_is_present   = .false.
    is_present             = .false.

    IF (present(cut_at)) THEN
       copy_cut = cut_at
    ELSE
       copy_cut = ' '
    ENDIF

    keyword_is_present = check_for(line,keyword,copy_cut)

    IF (present(and)) THEN
        second_word_is_present = check_for(line,and,copy_cut)
        keyword_is_present     = (keyword_is_present .and. second_word_is_present)
    ENDIF

    IF (present(alias) .and. .not. is_present) THEN
        alias_is_present = check_for(line,alias,copy_cut)
    ENDIF

    IF (present(but_not)) THEN
        exclusion_is_present = check_for(line,but_not,copy_cut)
    ENDIF

    is_present = (keyword_is_present .or. alias_is_present) .and. (.not. exclusion_is_present)


    CONTAINS

    ! ==--------------------------------------------------------------==
    FUNCTION check_for(line,keyword,cut_at) &
    RESULT (is_present)
 
      CHARACTER(len=*), INTENT(in)         :: line, keyword, cut_at
 
      CHARACTER(len=input_string_len)      :: next_words, cut_word
 
      LOGICAL                              :: is_present
      INTEGER                              :: first, last, cut_here

      is_present = .false.

      SELECT CASE(cut_at)
      CASE(" ")
         next_words = line
         full_keyword_loop: DO
            CALL xstring(next_words,first,last)
            IF (last - first < 0) THEN
               EXIT full_keyword_loop
            ELSEIF (trim(adjustl(next_words(first:last))) == trim(adjustl(keyword))) THEN
               is_present = .true.
               EXIT full_keyword_loop
            ENDIF
            next_words = next_words(last+1:)
         END DO full_keyword_loop
      CASE DEFAULT
         next_words = line
         cut_keyword_loop: DO
            CALL xstring(next_words,first,last)
            !
            ! Copy cut word
            !
            cut_word = trim(adjustl(next_words(first:last)))
            !
            ! Find index where to cut
            !
            cut_here = INDEX(cut_word,trim(adjustl(cut_at)))
            !
            ! Cut if necessary
            !
            IF (cut_here > 0) THEN
                cut_word = trim(adjustl(cut_word(1:cut_here-1)))
            ENDIF
            ! 
            IF (last - first < 0) THEN
               EXIT cut_keyword_loop
            ELSEIF (trim(adjustl(cut_word)) == trim(adjustl(keyword))) THEN
               is_present = .true.
               EXIT cut_keyword_loop
            ENDIF
            next_words = next_words(last+1:)
         END DO cut_keyword_loop
      END SELECT
 
    END FUNCTION check_for
    ! ==--------------------------------------------------------------==
  END FUNCTION keyword_contains
  ! ==================================================================
  FUNCTION index_of_delimiter(line,keyword,delimiter) &
  RESULT (delimiter_index)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), INTENT(in)         :: line, keyword, delimiter

    CHARACTER(len=input_string_len)      :: next_words, cut_word
 
    INTEGER                              :: first, last, delimiter_index, &
                                            running_index, next_index, cut_here

    delimiter_index = 0
    running_index   = 0

    !
    ! Loop to find the keyword
    !
    next_words = line
    find_keyword_loop: DO
       CALL xstring(next_words,first,last)
       running_index = running_index + first
       !
       ! Copy cut word
       !
       cut_word = trim(adjustl(next_words(first:last)))
       !
       ! Find index where to cut
       !
       cut_here = INDEX(cut_word,trim(adjustl(delimiter)))
       !
       ! Cut if necessary
       !
       IF (cut_here > 0) THEN
           cut_word = trim(adjustl(cut_word(1:cut_here-1)))
       ENDIF
       ! 
       IF (last - first < 0) THEN
          EXIT find_keyword_loop
       ELSEIF (trim(adjustl(cut_word)) == trim(adjustl(keyword))) THEN
          !
          ! Find out whether the delimiter is hidden in the keyword (keyword=)
          ! or whether we need to go on searching (keyword  = )
          !
          IF (cut_here > 0) THEN
             delimiter_index = cut_here + running_index
          ELSE
             !
             ! Take the next possible position of the delimiter
             !
             next_index = INDEX(next_words(first:),trim(adjustl(delimiter)))
             IF (next_index > 0) THEN
                delimiter_index = next_index + running_index
             ELSE
                !
                ! In case the delimiter cannot be found
                !
                delimiter_index = 0
             ENDIF
          ENDIF
          EXIT find_keyword_loop
       ENDIF
       next_words    = next_words(last+1:)
       running_index = running_index + (last - first)
    END DO find_keyword_loop
 
    ! ==--------------------------------------------------------------==
  END FUNCTION index_of_delimiter
  ! ==================================================================
  SUBROUTINE readsr(string,iin,iout,fptnum,error)
    ! ==--------------------------------------------------------------==
    ! == READ A FULL PRECISION (64 BIT) FLOATING POINT NUMBER         ==
    ! == FROM A STRING                                                ==
    ! == NOTE: ACCEPTED FORM IS                                       ==
    ! ==    (SIGN)___.___E(SIGN)__/(SIGN)___.___E(SIGN)__             ==
    ! ==    THE EXPONENT PORTION IS OPTIONAL AS IS THE DECIMAL POINT. ==
    ! ==    A FIELD STARTED BY E IS ASSUMED TO HAVE MANTISSA OF 1.    ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   STRING*(*) Character string                                ==
    ! ==   IIN  Index of the first read character in STRING           ==
    ! == OUTPUT:                                                      ==
    ! ==   IOUT Index of the last read character (STRING(IIN:IOUT))   ==
    ! ==   FPTNUM Full precision floating point number                ==
    ! ==   ERROR .TRUE. if error occured                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: string
    INTEGER                                  :: iin, iout
    REAL(real_8)                             :: fptnum
    LOGICAL                                  :: error

    REAL(real_8), PARAMETER                  :: one = 1._real_8, &
                                                ten = 10._real_8 , &
                                                zero = 0._real_8

    CHARACTER(len=80)                        :: iv
    INTEGER                                  :: i1, i2, Index, inv, isexp, j, &
                                                length, nexp
    LOGICAL                                  :: digit, exp, point, slash, &
                                                tsign
    REAL(real_8)                             :: scale, sign, v1

    CALL xstring(string(iin:LEN(string)),i1,i2)
    iout=iin+i2
    iv=string(iin+i1-1:iin+i2-1)
    length=i2-i1+1
    error=.FALSE.
    inv    = 0
    fptnum = zero
    IF (length.EQ.0) THEN
       error=.TRUE.
       RETURN
    ENDIF
    v1=zero
    sign=one
    scale=one
    digit = .FALSE.
    exp   = .FALSE.
    point = .FALSE.
    tsign = .FALSE.
    slash = .FALSE.
    isexp = 1
    nexp  = 0
    ! ==--------------------------------------------------------------==
    DO j=1,length
       ! 
       ! --- IDENTIFY THE J-TH CHARACTER ---
       ! 
       IF (iv(j:j).EQ.'0') THEN
          Index=0
       ELSEIF (iv(j:j).EQ.'1') THEN
          Index=1
       ELSEIF (iv(j:j).EQ.'2') THEN
          Index=2
       ELSEIF (iv(j:j).EQ.'3') THEN
          Index=3
       ELSEIF (iv(j:j).EQ.'4') THEN
          Index=4
       ELSEIF (iv(j:j).EQ.'5') THEN
          Index=5
       ELSEIF (iv(j:j).EQ.'6') THEN
          Index=6
       ELSEIF (iv(j:j).EQ.'7') THEN
          Index=7
       ELSEIF (iv(j:j).EQ.'8') THEN
          Index=8
       ELSEIF (iv(j:j).EQ.'9') THEN
          Index=9
       ELSEIF (iv(j:j).EQ.'/') THEN
          Index=10
       ELSEIF (iv(j:j).EQ.'.') THEN
          Index=11
       ELSEIF (iv(j:j).EQ.'+') THEN
          Index=12
       ELSEIF (iv(j:j).EQ.'-') THEN
          Index=13
       ELSEIF (iv(j:j).EQ.'E') THEN
          Index=14
       ELSEIF (iv(j:j).EQ.'e') THEN
          Index=14
       ELSEIF (iv(j:j).EQ.'D') THEN
          Index=15
       ELSEIF (iv(j:j).EQ.'d') THEN
          Index=15
       ELSE
          ! ILLEGAL CHARACTER FOUND
          error=.TRUE.
          RETURN
       ENDIF
       ! 
       ! --- TEST IF THIS CHARACTER IS A NON-NUMERAL ---
       ! 
       IF (Index.LE.9) THEN
          ! FIND A DIGIT
          digit=.TRUE.
          tsign=.FALSE.
          ! TEST IF PART OF EXPONENT
          IF (exp) THEN
             ! ADD DIGIT TO EXPONENT
             nexp=10*nexp+Index
             CYCLE
          ENDIF
          ! ADD DIGIT TO MANTISSA
          inv=inv+1
          IF (point) THEN
             scale=scale/ten
             fptnum=fptnum+REAL(Index,kind=real_8)*scale
          ELSE
             fptnum=ten*fptnum+REAL(Index,kind=real_8)
          ENDIF
          CYCLE
       ENDIF
       ! 
       ! --- PROCESS NON-NUMERALS CHARACTERS ---
       IF (Index.EQ.10) THEN
          GOTO 120
       ELSEIF (Index.EQ.11) THEN
          GOTO 125
       ELSEIF (Index.EQ.12) THEN
          GOTO 130
       ELSEIF (Index.EQ.13) THEN
          GOTO 135
       ELSEIF (Index.EQ.14.OR.index.EQ.15) THEN
          GOTO 145
       ENDIF
       ! 
       ! --- SLASH DETECTED (FIELD EXPRESSED AS A FRACTION) NUMERATOR COMP
       ! 
120    CONTINUE
       IF (slash) THEN
          error=.TRUE.
          RETURN
       ENDIF
       slash=.TRUE.
       v1=fptnum*sign*(ten**(isexp*nexp))
       fptnum=zero
       tsign=.FALSE.
       digit=.FALSE.
       inv=0
       IF (v1.EQ.zero) go to 155
       sign=one
       isexp=1
       exp=.FALSE.
       scale=one
       point=.FALSE.
       nexp=0
       CYCLE
       ! 
       ! --- DECIMAL POINT DETECTED ---
       ! 
125    CONTINUE
       IF (exp) THEN
          error=.TRUE.
          RETURN
       ENDIF
       point=.TRUE.
       CYCLE
       ! 
       ! --- PLUS TSIGN DETECTED, TEST IF START OF MANTISSA OR EXPONENT ---
       ! 
130    CONTINUE
       IF (tsign) THEN
          error=.TRUE.
          RETURN
       ELSE
          tsign=.TRUE.
       ENDIF
       IF (exp) THEN
          IF (nexp.NE.0) THEN
             error=.TRUE.
             RETURN
          ENDIF
          isexp=1
       ELSE
          IF (inv.NE.0) THEN
             error=.TRUE.
             RETURN
          ENDIF
          sign=one
       ENDIF
       CYCLE
       ! 
       ! --- MINUS TSIGN, TEST IF START OF MANTISSA OF EXPONENT ---
       ! 
135    CONTINUE
       IF (tsign) THEN
          error=.TRUE.
          RETURN
       ELSE
          tsign=.TRUE.
       ENDIF
       IF (exp) THEN
          IF (nexp.NE.0) THEN
             error=.TRUE.
             RETURN
          ENDIF
          isexp=-1
       ELSE
          IF (inv.NE.0) THEN
             error=.TRUE.
             RETURN
          ENDIF
          sign=-one
       ENDIF
       CYCLE
       ! 
       ! --- E, D OR EMBEDDED + OR - STARTS EXPONENT FIELD ---
       ! 
145    CONTINUE
       IF (exp) THEN
          error=.TRUE.
          RETURN
       ENDIF
       exp=.TRUE.
       ! A FIELD STARTED BY E IS ASSUMED TO HAVE MANTISSA OF 1.
       IF (inv.EQ.0) fptnum=one
       tsign=.FALSE.
       CYCLE
    END DO
    ! ==--------------------------------------------------------------==
    IF (digit) THEN
       fptnum=fptnum*sign*(ten**(isexp*nexp))
    ELSE
       ! NO DIGIT DETECTED
       error=.TRUE.
       RETURN
    ENDIF
    ! 
    ! --- THE NUMERATOR IS FINISHED, TEST FOR NO DENOMINATOR ---
    ! 
    IF (v1.EQ.zero) go to 155
    IF (fptnum.EQ.zero) fptnum=one
    fptnum=v1/fptnum
155 CONTINUE
    IF (digit) THEN
       v1=fptnum
    ELSE
       ! NO DIGIT DETECTED
       error=.TRUE.
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE readsr
  ! ==================================================================
  SUBROUTINE readsi(string,iin,iout,intnum,error)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: string
    INTEGER                                  :: iin, iout, intnum
    LOGICAL                                  :: error

    REAL(real_8)                             :: fptnum

    CALL readsr(string,iin,iout,fptnum,error)
    intnum=NINT(fptnum)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE readsi
  ! ==================================================================
  SUBROUTINE xstring(string,ia,ie)
    ! ==--------------------------------------------------------------==
    ! == Give index of the first and last non-blank character         ==
    ! == blank characters are " " and CHAR(9) (tab)                   ==
    ! == more blank characters are char(13) (cr) and CHAR(10) (newl)  ==
    ! == STRING (INPUT)                                               ==
    ! == IA,IE  (OUTPUT)                                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: string
    INTEGER                                  :: ia, ie

    INTEGER                                  :: i, slen

    slen=LEN(string)
    ia=1
    DO i=1,slen
       IF (string(i:i).NE.' '.AND.&
            string(i:i).NE.CHAR(9).AND.      & ! tab character
            string(i:i).NE.CHAR(10).AND.     & ! newline character
            string(i:i).NE.CHAR(13) ) THEN   ! carriage return character
          ia=i
          GOTO 10
       ENDIF
    ENDDO
10  CONTINUE
    DO i=ia,slen
       IF (string(i:i).EQ.' '.OR.&
            string(i:i).EQ.CHAR(0).OR. & ! \0 character
            string(i:i).EQ.CHAR(9)) THEN ! tab character
          ie=i-1
          GOTO 20
       ENDIF
    ENDDO
    ie=slen
20  CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xstring
  ! ==================================================================
END MODULE readsr_utils
