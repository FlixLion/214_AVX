
default rel
global	as2v_complex
global	_as2v_complex


;addsig2vol_2(zeros(3000,2),single(zeros(3,1)),single(zeros(3,1)),single(zeros(3,1)),single(zeros(1,1),single(zeros(1,1),single(zeros(3,1),uint32([10 10 10]),zeros([10 10 10]))  );
;addsig2vol_2(repmat([1:3000]',[1 2]),single(ones(3,1)),single(ones(3,2)),single(ones(3,1)),single(ones(1,1)),single(ones(1,1)),single(ones(3,1)),uint32([10 10 10]),zeros([10 10 10])  );

[segment .data align=16]
;;;;variables
	WIDTH_BUFFER dd 0
	MASK_BUFFER dd 4294967295
	RSP_SAVE dq 0

	
[segment .bss align=16]
 XMM6_SAVE resq 2 ; in MASM: 2 dup (?)
 XMM7_SAVE resq 2 ; in MASM: 2 dup (?)


[segment .text align=16]

;void as2v(double* bild, double* AScan,int n_AScan, double* buffer, double* index_vol, int n_index, double* rec_pos,double* em2pix, double *speed, double* resolution,double* time intervall,double* AScan_complex, double* buffer_complex, double* bild)+68

;;;;;;
;old ;                  new			        			   newest 64                        newest 64
;;;;;;		            	;;;;;				         	   	;;;;;; (reording to get around struct align)
;28=*bild								                            	rcx+0   *bild	                 rcx+0
;32=*Ascan								                           	rcx+8   *bild_compl                rcx+8
;36=n_Ascan				        	   	        ->*bild_compl	rcx+16  *Ascan
;40=*buffer								                          	rcx+24  *Ascan_complex
;44=*pixel_index      ->  *pixel_vect[1 1 1]	    		rcx+32  *buffer	
;48=n_index           ->  n_pixX		   -> *buff_cmpl	rcx+40  *buffer_compl
;52=*receiverpos							                      	rcx+48  *IMAGE_SUM
;56=*dist_em_pixel    ->  *senderpos			          	rcx+56  *IMAGE_SUM_compl
;60=*speed								                           	rcx+64
;64=*resolut							                          	rcx+72
;68=*time_int						                          		rcx+80
;72=*Ascan_complex				                  		    	rcx+88
;76=*buff_complex					              -> *imag_sum	rcx+96
;80=*bild_compl					              -> *imag_sum_c	rcx+104
;84=                  ->  n_pixY				            	rcx+112
;88=                  ->  n_pixZ			               	rcx+116
;92=                  ->  *IMAGE_SUM	      -> nAscan	rcx+120
;96=                  ->  *IMAGE_SUM_compl  -> n_pixX	rcx+124
;100=                 ->  quadwordbuffer    -> NULL 	rcx+128
	
    

as2v_complex:	
_as2v_complex:

	
	push rbp
	push rsi
	push rdi
	push rbx	 ;+12 (diff:esp+8)

	;push eax	 ;buffer for width  +8 (diff esp+12)
	;push eax	 ;buffer for CTW    +4 (diff esp+16)
 

	;save rsp
	mov [rel RSP_SAVE],qword rsp
	
	
	;save XMM6-7 (callee convention WIN64(?))
	movups [rel XMM6_SAVE],xmm6
	movups [rel XMM7_SAVE],xmm7

	;transfer via rcx because first parameter on MS-64 convention and 4. on linux,so pass the pointer as 1. & 4. 
	mov rsp,rcx

 	fninit		 ;init fpu


   ;calc width fpu
	;FNSTCW [esp+4]   ;buffering in stack (16bit!) + FWAIT to ensure save is complete (FSTCW without wait)
	;mov ax,[esp+4]
	;and ax,0F3FFh    ;delete 2 bits
	;or  ax,0000100000000000b ;round up +inf
	;mov [esp+6],ax   ;new CTRWORD
	;FLDCW [esp+6]

	;fldz                 ;load zero
	;fld1                 ;load 1
	;mov eax,[rsp+80]     ;time interval
	    ;fld qword [eax]     ;DOUBLE = qword
	;fld dword [eax]      ;SINGLE = dword
	;mov eax,[rsp+72]     ;resolution
	     ;fld qword [eax]       ;DOUBLE = qword
	;fld dword [eax]      ;SINGLE = dword
	;mov eax,[rsp+64]     ;speed
		;fld qword [eax]       ;DOUBLE = qword
	;fld dword [eax]      ;SINGLE = dword

	;fdivp st1,st0        ; res/speed; destination/source ,pop
	;fdivrp st1,st0       ; ergeb/timeinterv; source/dest ,pop -> width now in st0

	;fcom st0,st1         ; z.b.: 0.4 cmp 1  -> lower
	;fcmovb st0,st2       ; transfer st2 to st0 if below -> width<1 -> interpol not sum!

	;mov [esp+8],dword 2
	;fidiv dword [esp+8]  ; halbieren
	;fistp dword [esp+8]  ; width to int32 (ceiling round)

	;FLDCW [esp+4]        ;restore old CTRWORD
	;fstp    st0
	;fstp    st0          ; -> fpu stack emtpy now
       ;end calc width fpu
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     ;;;calc width with SSE1
      mov eax,1
      cvtsi2ss xmm0,eax
      mov eax,2
      cvtsi2ss xmm1,eax
      mov rax,[rsp+80]	   ;*time interval
      movss xmm2,[rax]
      mov rax,[rsp+72]	   ;*resolution
      movss xmm3,[rax]
      mov rax,[rsp+64]	   ;*speed
      movss xmm4,[rax]

      divss xmm3,xmm4	   ;res/speed -> time
      divss xmm3,xmm2	   ;time/timeinterv. -> number of samples -> width
      addss xmm3,xmm0	   ;width+1  (against floor)
      divss xmm3,xmm1	   ;width/2
      cvttss2si eax,xmm3   ;to int32 (floor ROUNDING)
     ;;;end  calc width with SSE1

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;;width to big?
	mov ecx,[rsp+120]  ;n_Ascan
	shr ecx,1	  ;n_aScan/2
	sub ecx,2
	cmp eax,ecx
	jl xs_int_check   ;width < (n_ascan/2)-2
	mov eax,ecx	  ;set Width to maximum ((n_ascan/2)-2)

      ;xsum or interpol?
      xs_int_check:
	cmp eax,0	  ;width 0 check (and negative values!)
	jg xsum
	jmp interpol
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;;Xsum path
	mov dword [rel WIDTH_BUFFER],eax   ;width buffer

     xsum:
	mov rbx,[rsp+24]  ;*Buffer
	mov rdx,[rsp+8]  ;*AScan
	mov edi,eax	      ;edi=width
	mov rsi,[rsp+120]  ;n_AScan
			  ;eax = width

	fld qword [rdx]   ;AScan[0]
      b0:
	add rdx,8
	fadd qword [rdx]
	dec eax     	  ;eax = width
	jnz b0		  ;AScan[width]

	mov ecx,edi	      ;Width
	fst qword [rbx]   ;buffer[0]

      b1:
	add rdx,8
	fadd qword [rdx]
	add rbx,8
	fst qword [rbx]   ;buffer[1]
	dec ecx
	jnz b1

	mov rax,rdx	  ; edx ->vorauseilend um width
	xor rcx,rcx ; clear RCX for later mixed 32-64 bit add
	mov ecx,edi	  ; width
	shl ecx,4	  ;imul ecx,16       ; 2*width ->addr
	add ecx,8
	sub rax,rcx	  ;hinterher um 2*width +1

	mov ecx,esi	  ;n_AScan
	sub ecx,edi	  ; -width
	sub ecx,edi	  ; -width
	sub ecx,1	  ;n_AScan-2*Width-1

	mov rbp,8
      b2:
	fadd qword [rdx+rbp] ;vorauseilend
	fsub qword [rax+rbp] ;hinterher
	fst qword [rbx+rbp]  ;buffer[ebx]
	add rbp,8  
	sub ecx,1
	jnz b2		     ; n_AScan - width

	add rdx,rbp
	add rax,rbp
	add rbx,rbp

	mov ecx,edi   ;width
      b3:
       fsub qword [rax]
	add rax,8	  ;hinterher
       fst qword [rbx]	  ;buffer[rbx]
	add rbx,8

	dec ecx
	jnz b3		  ;n_AScan

	fstp  st0	  ;pop st0
     ;end create xsum buffer_real

     ;l_interpol_real (range 5x)
	mov rdi,[rsp+24]  ;*Buffer_read[0]    ;;;;;[rsp+8] *Ascan[0]
	xor rcx,rcx ; clear RCX for later mixed 32-64 bit add
	mov ecx,[rsp+120]  ;n_AScan
	dec ecx
	shl ecx,3	  ;double ptr *ascan[2999]
	add rdi,rcx	  ;double ptr *ascan[2999]

	mov rsi,[rsp+24]  ;*Buffer_write[0]
	mov rbp,rsi	  ;*Buffer_write[0]
	mov rdx,[rsp+120]  ;n_AScan
	xor rax,rax ; clear RAX for later mixed 32-64 bit add
	mov eax,5	  ;interp_range!!!
	mul edx 	  ;eax*edx  result 64bit: [EDX EAX] !!!!
	dec eax
	shl eax,3	  ;double ptr *buffer[14999]
	add rsi,rax	  ;double ptr *buffer[14999]

	;prefill FPU Stack
	fld1
	fld1
	fadd st1,st0
	fldz
	fldz
	fldz
	fadd st0,st4  ;st0=2
	fadd st0,st4  ;st0=4
	fadd st0,st0  ;st0=8
	fadd st0,st4  ;st4=2, st3=1; st2=0; st1=0; st0=10

	fdiv st4,st0	   ;0,2 = st4; st3 = ???(buffer), st2 = ???(buffer)  st1 = ???(buffer), st0 = ???(buffer)

	;last 4pixel position
	fld qword [rdi]    ;*Ascan[2999]
	fmul st0,st5	   ;0,2 * (*Ascan[2999])
	fst qword [rsi]    ; *Buffer[14999]  0.2
	fst st1
	fadd st0,st1
	fst qword [rsi-8]  ; *Buffer[14998]  0.4
	fadd st0,st1
	fst qword [rsi-16] ; *Buffer[14997]  0.6
	fadd st0,st1	   ; 0,2 * (*Ascan[2999]) + 0,2 * (*Ascan[2999])
	fstp qword [rsi-24]; *Buffer[14996]  0.8

	sub rsi,32

       ;FILL loop
       loop_interp_real:

	fld qword [rdi-8] ;*Ascan[2998]
	fld qword [rdi]   ;*Ascan[2999]
	fst qword [rsi]   ;*Buffer[14995]

	fmul st0,st6	  ; *Ascan[2999] *0.2
	fxch st1,st0
	fmul st0,st6	  ; *Ascan[2998] *0.2 -> st0=0.2*2998 st1=0.2*2999 st2=? st3=? st4=? st5=? st6=0.2

	fst st2
	fst st3
	fst st4
	fst st5
	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st2,st0
	fadd st3,st0
	fadd st2,st0

	fxch st0,st1

	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st4,st0
	fadd st5,st0
	faddp st5,st0

	fstp st0	    ; st0=1 st1=2 st2=3 st3=4 st4=0.2
	fst qword [rsi-32]
	fxch st0,st1
	fst qword [rsi-24]
	fxch st0,st2
	fst qword [rsi-16]
	fxch st0,st3
	fst qword [rsi-8]

	sub rsi,40
	sub rdi,8

	cmp rsi,rbp	    ;exit because *Buffer[0] reached
	jne loop_interp_real

	;beginnig range
	fld qword [rdi]     ;*Ascan[0]
	fstp qword [rsi]    ;*Buffer[0]

	fstp st0
	fstp st0
	fstp st0
	fstp st0
	fstp st0	   ; FPU stack empty
	;;;;;;;;;;;;;;;;;;;Buffer real filled

      ;ascan complex?
      mov rax,[rsp+88]	  ;*AScan_complex
      cmp rax,0
      jne buff_complex1
      jmp fill_real_bild  ; NO Complex AScan

     ;create buffer complex
     buff_complex1:
	mov ecx,[rel WIDTH_BUFFER]   ;Width
	mov rbx,[rsp+40]  ;*Buffer_complex
	mov rdx,[rsp+88]  ;*AScan_complex

	fld qword [rdx]   ;AScan_complex[0]
      b0_c:
	add rdx,8
	fadd qword [rdx]
	dec ecx
	jnz b0_c	  ;AScan[width]

	mov ecx,[rel WIDTH_BUFFER]   ;Width
	fst qword [rbx]   ;buffer[0]

      b1_c:
	add rdx,8
	fadd qword [rdx]
	add rbx,8
	fst qword [rbx]   ;buffer[1]
	dec ecx
	jnz b1_c


	mov rax,rdx	  ; rdx ->vorauseilend um width
	xor rcx,rcx ; clear RCX for later mixed 32-64 bit add
	mov ecx,[rel WIDTH_BUFFER]   ; width
	shl ecx,4	  ;imul ecx,16       ; 2*width ->addr
	add ecx,8
	sub rax,rcx	  ;hinterher um 2*width +1

	mov ecx,[rsp+120]  ;n_AScan
	sub ecx,[rel WIDTH_BUFFER]   ; -width
	sub ecx,[rel WIDTH_BUFFER]   ; -width
	sub ecx,1	  ;n_AScan-2*Width-1

	mov rbp,8
      b2_c:
	fadd qword [rdx+rbp] ;vorauseilend
	fsub qword [rax+rbp] ;hinterher
	fst qword [rbx+rbp]  ;buffer[ebx]
	add rbp,8
	sub ecx,1
	jnz b2_c	     ; n_AScan - width

	add rdx,rbp
	add rax,rbp
	add rbx,rbp

	mov ecx,[rel WIDTH_BUFFER]   ;width
      b3_c:
	add rax,8	  ;hinterher
	fsub qword [rax]
	add rbx,8
	fst qword [rbx]   ;buffer[rbx]
	dec ecx
	jnz b3_c	  ;n_AScan

	fstp  st0	  ;pop st0
     ;end create xsum buffer_complex

     ;l_interpol_COMPLEX (range 5x)
  xor rcx,rcx ; clear RCX for later mixed 32-64 bit add
	mov rdi,[rsp+40]  ;*Buffer_read[0]    ;;;;;[rsp+8] *Ascan[0]
	mov ecx,[rsp+120]  ;n_AScan
	dec ecx
	shl ecx,3	  ;double ptr *ascan[2999]
	add rdi,rcx	  ;double ptr *ascan[2999]

	mov rsi,[rsp+40]  ;*Buffer_write[0] as countdown variable
	mov rbp,rsi	      ;*Buffer_write[0] to check later against

	xor rax,rax ; clear RAX for later mixed 32-64 bit add
	mov edx,[rsp+120]  ;n_AScan
	mov eax,5	  ;interp_range!!!
	mul edx 	  ;eax*edx  result 64bit: [EDX EAX] !!!!
	dec eax
	shl eax,3	  ;double ptr *buffer[14999]
	add rsi,rax	  ;double ptr *buffer[14999]

	;prefill FPU Stack
	fld1
	fld1
	fadd st1,st0
	fldz
	fldz
	fldz
	fadd st0,st4  ;st0=2
	fadd st0,st4  ;st0=4
	fadd st0,st0  ;st0=8
	fadd st0,st4  ;st4=2, st3=1; st2=0; st1=0; st0=10

	fdiv st4,st0	 ;0,2 = st4; st3 = ???(buffer), st2 = ???(buffer)  st1 = ???(buffer), st0 = ???(buffer)

	;last 4pixel position
	fld qword [rdi]    ;*Ascan[2999]
	fmul st0,st5	   ;0,2 * (*Ascan[2999])
	fst qword [rsi]    ; *Buffer[14999]  0.2
	fst st1
	fadd st0,st1
	fst qword [rsi-8]  ; *Buffer[14998]  0.4
	fadd st0,st1
	fst qword [rsi-16] ; *Buffer[14997]  0.6
	fadd st0,st1	   ; 0,2 * (*Ascan[2999]) + 0,2 * (*Ascan[2999])
	fstp qword [rsi-24]; *Buffer[14996]  0.8

	sub rsi,32

       ;FILL loop
       loop_interp_complex:

	fld qword [rdi-8] ;*Ascan[2998]
	fld qword [rdi]   ;*Ascan[2999]
	fst qword [rsi]   ;*Buffer[14997]

	fmul st0,st6	  ; *Ascan[2999] *0.2
	fxch st1,st0
	fmul st0,st6	  ; *Ascan[2998] *0.2 -> st0=0.2*2998 st1=0.2*2999 st2=? st3=? st4=? st5=? st6=0.2

	fst st2
	fst st3
	fst st4
	fst st5
	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st2,st0
	fadd st3,st0
	fadd st2,st0

	fxch st0,st1

	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st4,st0
	fadd st5,st0
	faddp st5,st0

	fstp st0	; st0=1 st1=2 st2=3 st3=4 st4=0.2
	fst qword [rsi-32]
	fxch st0,st1
	fst qword [rsi-24]
	fxch st0,st2
	fst qword [rsi-16]
	fxch st0,st3
	fst qword [rsi-8]

	sub rsi,40
	sub rdi,8

	cmp rsi,rbp	;exit because *Buffer[0] reached
	jne loop_interp_complex

	;beginnig range
	fld qword [rdi]     ;*Ascan[0]
	fstp qword [rsi]    ;*Buffer[0]

	fstp st0
	fstp st0
	fstp st0
	fstp st0
	fstp st0	   ; FPU stack empty
	;;;;;;;;;;;;;;;;;;;Buffer COMPLEX fill

     jmp fill_real_bild
    ;;;;;end xsum preparation


    ;;;;;;;;;;;;;;;;;;;begin interpol
    interpol:

      ;l_interpol_real (range 5x)
	mov rdi,[rsp+8]  ;*Ascan[0]
  xor rcx,rcx ; clear RCX for later mixed 32-64 bit add
	mov ecx,[rsp+120]  ;n_AScan
	dec ecx
	shl ecx,3	  ;double ptr *ascan[2999]
	add rdi,rcx	  ;double ptr *ascan[2999]

	mov rsi,[rsp+24]  ;*Buffer_write[0] as countdown variable
	mov rbp,rsi	  ;*Buffer_write[0] to check against on end

	mov edx,[rsp+120]  ;n_AScan
	or rax,rax ; clear RAX for later mixed 32-64 bit add
	mov eax,5	  ;interp_range!!!
	mul edx 	  ;eax*edx  result 64bit: [EDX EAX] !!!!
	dec eax
	shl eax,3	  ;double ptr *buffer[14999]
	add rsi,rax	  ;double ptr *buffer[14999]

	;prefill FPU Stack
	fld1
	fld1
	fadd st1,st0
	fldz
	fldz
	fldz
	fadd st0,st4  ;st0=2
	fadd st0,st4  ;st0=4
	fadd st0,st0  ;st0=8
	fadd st0,st4  ;st4=2, st3=1; st2=0; st1=0; st0=10

	fdiv st4,st0	 ;0,2 = st4; st3 = ???(buffer), st2 = ???(buffer)  st1 = ???(buffer), st0 = ???(buffer)

	;last 4pixel position
	fld qword [rdi]    ;*Ascan[2999]
	fmul st0,st5	   ;0,2 * (*Ascan[2999])
	fst qword [rsi]    ; *Buffer[14999]  0.2
	fst st1
	fadd st0,st1
	fst qword [rsi-8]  ; *Buffer[14998]  0.4
	fadd st0,st1
	fst qword [rsi-16] ; *Buffer[14997]  0.6
	fadd st0,st1	   ; 0,2 * (*Ascan[2999]) + 0,2 * (*Ascan[2999])
	fstp qword [rsi-24]; *Buffer[14996]  0.8

	sub rsi,32

       ;FILL loop
       loop_interp_only_real:

	fld qword [rdi-8] ;*Ascan[2998]
	fld qword [rdi]   ;*Ascan[2999]
	fst qword [rsi]   ;*Buffer[14997]

	fmul st0,st6	  ; *Ascan[2999] *0.2
	fxch st1,st0
	fmul st0,st6	  ; *Ascan[2998] *0.2 -> st0=0.2*2998 st1=0.2*2999 st2=? st3=? st4=? st5=? st6=0.2

	fst st2
	fst st3
	fst st4
	fst st5
	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st2,st0
	fadd st3,st0
	fadd st2,st0

	fxch st0,st1

	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st4,st0
	fadd st5,st0
	faddp st5,st0

	fstp st0	; st0=1 st1=2 st2=3 st3=4 st4=0.2
	fst qword [rsi-32]
	fxch st0,st1
	fst qword [rsi-24]
	fxch st0,st2
	fst qword [rsi-16]
	fxch st0,st3
	fst qword [rsi-8]

	sub rsi,40
	sub rdi,8

	cmp rsi,rbp	;exit because *Buffer[0] reached
	jne loop_interp_only_real

	;beginning range
	fld qword [rdi]     ;*Ascan[0]
	fstp qword [rsi]    ;*Buffer[0]

	fstp st0
	fstp st0
	fstp st0
	fstp st0
	fstp st0	   ; FPU stack empty
	;;;;;;;;;;;;;;;;;;;Buffer real filled

      ;ascan complex?
      mov rax,[rsp+88]	  ;*AScan_complex
      cmp rax,0
      jne buff_complex2
      jmp fill_real_bild  ; NO Complex AScan

     ;l_interpol_COMPLEX (range 5x)
     buff_complex2:
	mov rdi,[rsp+88]  ;*Ascan_complex_read[0]
	xor rcx,rcx ; clear RCX for later mixed 32-64 bit add
	mov ecx,[rsp+120]  ;n_AScan
	dec ecx
	shl ecx,3	  ;double ptr *ascan[2999]
	add rdi,rcx	  ;double ptr *ascan[2999]

	mov rsi,[rsp+40]  ;*Buffer_write[0] count down variable
	mov rbp,rsi	  ;*Buffer_write[0] to check against
  
  xor rax,rax ; clear rax for later mixed 32-64 bit add
	mov edx,[rsp+120]  ;n_AScan
	mov eax,5	  ;interp_range!!!
	mul edx 	  ;eax*edx  result 64bit: [EDX EAX] !!!!
	dec eax
	shl eax,3	  ;double ptr *buffer[14999]
	add rsi,rax	  ;double ptr *buffer[14999]

	;prefill FPU Stack
	fld1
	fld1
	fadd st1,st0
	fldz
	fldz
	fldz
	fadd st0,st4  ;st0=2
	fadd st0,st4  ;st0=4
	fadd st0,st0  ;st0=8
	fadd st0,st4  ;st4=2, st3=1; st2=0; st1=0; st0=10

	fdiv st4,st0	 ;0,2 = st4; st3 = ???(buffer), st2 = ???(buffer)  st1 = ???(buffer), st0 = ???(buffer)

       ;last 4pixel position
	fld qword [rdi]    ;*Ascan[2999]
	fmul st0,st5	   ;0,2 * (*Ascan[2999])
	fst qword [rsi]    ; *Buffer[14999]  0.2
	fst st1
	fadd st0,st1
	fst qword [rsi-8]  ; *Buffer[14998]  0.4
	fadd st0,st1
	fst qword [rsi-16] ; *Buffer[14997]  0.6
	fadd st0,st1	   ; 0,2 * (*Ascan[2999]) + 0,2 * (*Ascan[2999])
	fstp qword [rsi-24]; *Buffer[14996]  0.8

	sub rsi,32

       ;FILL loop
       loop_interp_only_complex:

	fld qword [rdi-8] ;*Ascan[2998]
	fld qword [rdi]   ;*Ascan[2999]
	fst qword [rsi]   ;*Buffer[14997]

	fmul st0,st6	  ; *Ascan[2999] *0.2
	fxch st1,st0
	fmul st0,st6	  ; *Ascan[2998] *0.2 -> st0=0.2*2998 st1=0.2*2999 st2=? st3=? st4=? st5=? st6=0.2

	fst st2
	fst st3
	fst st4
	fst st5
	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st2,st0
	fadd st3,st0
	fadd st2,st0

	fxch st0,st1

	fadd st2,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st3,st0
	fadd st4,st0
	fadd st5,st0
	fadd st4,st0
	fadd st5,st0
	faddp st5,st0

	fstp st0	; st0=1 st1=2 st2=3 st3=4 st4=0.2
	fst qword [rsi-32]
	fxch st0,st1
	fst qword [rsi-24]
	fxch st0,st2
	fst qword [rsi-16]
	fxch st0,st3
	fst qword [rsi-8]

	sub rsi,40
	sub rdi,8

	cmp rsi,rbp	;exit because *Buffer[0] reached
	jne  loop_interp_only_complex

	;beginning range
	fld qword [rdi]     ;*Ascan[0]
	fstp qword [rsi]    ;*Buffer[0]

	fstp st0
	fstp st0
	fstp st0
	fstp st0
	fstp st0	   ; FPU stack empty
	;;;;;;;;;;;;;;;;;;;Buffer COMPLEX fill


    jmp fill_real_bild
    ;;;;;;;;;;;;;;end interpol


 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;fill images
    ;;;;;;;;;;;;;;;;;;;;;;;;;
    fill_real_bild:

    ;init mm7 with zero
    ;mov eax,0           ;
    ;cvtsi2ss xmm7,eax   ; xmm7 = 0
    ;shufps   xmm7,xmm7,0; 0,0,0,0
    ;CVTPS2PI mm7,xmm7   ; mm7  = 0

    ;pixel_pos
    mov rdx,[rsp+32]	 ;*pixelpos
    movss  xmm7,[rdx+8]  ;pixelpos      0,0,0,Z
    shufps xmm7,xmm7,0	 ;pixelpos      Z,Z,Z,Z
    movlps xmm7,[rdx]	 ;pixelpos      Z(TRASH),Z,Y,X

    ;rec_pos
    mov rdx,[rsp+48]	 ;*receiverpos
    movss  xmm6,[rdx+8]  ;receiverpos   0,0,0,Z
    shufps xmm6,xmm6,0	 ;receiverpos   Z,Z,Z,Z
    movlps xmm6,[rdx]	 ;receiverpos   Z(TRASH),Z,Y,X

    ;sender_pos
    mov rdx,[rsp+56]	 ;*senderpos
    movss  xmm5,[rdx+8]  ;senderpos     0,0,0,Z
    shufps xmm5,xmm5,0	 ;senderpos     Z,Z,Z,Z
    movlps xmm5,[rdx]	 ;senderpos     Z(TRASH),Z,Y,X

    ;resolut
    mov rdx,[rsp+72]	 ;*resolut
    movss xmm4,[rdx]	 ;t,t,t,resolut

    ;factor = 1/(speed*timeinterval)
    mov eax,5		 ; 1 for interp_ratio 1 -> 5 for actual !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cvtsi2ss xmm3,eax	 ;xmm3=1
    mov rdx,[rsp+64]	 ;*speed
    movss xmm2,[rdx]
    mov rdx,[rsp+80]	 ;*timeinterval
    movss xmm1,[rdx]
    mulss xmm2,xmm1
    divss xmm3,xmm2	 ;xmm3 = 1/(speed*timeinterval)
    unpcklps xmm4,xmm3	 ; t,t,factor,resolut

    ;interleave in xmm4
    mov eax,-1		 ;4294967295
    cvtsi2ss xmm3,eax	 ;xmm3 = 1

    ;fixed interpol factor
    mov eax,5
    cvtsi2ss xmm1,eax	 ;xmm1 = 5

    mov eax,[rsp+120]	 ;n_Ascan
    cvtsi2ss xmm2,eax	 ;xmm2 = 3000

    mulss xmm2,xmm1	 ;3000*5 ->15000
    cvtss2si eax,xmm2;
    mov dword [rel WIDTH_BUFFER],eax ;15000
    unpcklps xmm3,xmm2	 ;xmm3 = t,t,15000,-1
    shufps xmm4,xmm3,68  ;0100 0100b -> xmm4 = 15000,-1,factor,resolut
    
    ;;;create complete factor buffer
    movaps xmm8,xmm4
	  shufps xmm8,xmm8,85  ;    0101 0101b -> factor,factor,factor,factor
    
    ;;; create 2*resolut buffer
    movss  xmm10,xmm4    ;    t,t,t,resolut
    shufps xmm10,xmm10,0 ;    0000 0000b -> resolut,resolut,reolut,resolut 
    addps xmm10,xmm10    ;    2*resolut,2*resolut,2*resolut,2*resolut  
    
    ;;fill up counters
    mov edi,[rsp+116]	 ; n_pixZ
    mov edx,[rsp+112]	 ; n_pixY
    mov ecx,[rsp+124]	 ; n_pixX
   
    ;;;fillup pointer to data
    mov rsi,[rsp+88]	 ; *ascan_complex (if not existent = NULL = 0)
 
    mov r10,[rsp+32]  ;*pixelpos
    mov r11,[rsp+24]	;*buffer_real[index_aktual]               NP
	  mov r12,[rsp+96]	;*image_sum_real
    mov r13,[rsp+0]	  ;*bild_real
     
    movss  xmm15,[rel MASK_BUFFER];	 xmm15 = ffffffffh
    shufps xmm15,xmm15,0   ;0 -> xmm15 = 4x ffffffffh
    
    
    ;;delta Counter
    mov rbp,0		 ; delta IMAGE

   Z_loop:

     Y_loop_init:
     mov edx,[rsp+112]	  ; n_pixY
     movlps xmm7,[r10]	  ;pixelpos     (t),Z(akt),Y(0),X(0_but_is_notimportant)

     Y_loop:
       ;pix z,y
       movaps xmm3,xmm7     ;          (trash),z_p,y_p,x_p
       shufps xmm3,xmm3,165 ; 10100101b -> z_p,z_p,y_p,y_p
       ;receiver z,y
       movaps xmm2,xmm6     ;     (trash),z_r,y_r,x_r
       shufps xmm2,xmm2,153 ; 10011001b -> (t),(t),z_r,y_r       psrldq xmm2,4        ;      SSE2 shift right by BYTES not BITS! t,t,z,y
       ;sender  z,y
       movaps xmm1,xmm5     ;     (trash),z_s,y_s,x_s
       shufps xmm1,xmm1,153 ; 10011001b -> (t),(t),z_s,y_s       psrldq xmm1,4        ;      SSE2 shift right by BYTES not BITS! t,t,z,y
       ;interleave
       unpcklps xmm2,xmm1   ;     z_s,z_r,y_s,y_r
    
       subps xmm3,xmm2	    ;     pix-senderpos pix-receiverpos
       mulps xmm3,xmm3	    ;     quadrieren

       movaps xmm2,xmm3     ;     z_s,z_r,y_s,y_r
       shufps xmm2,xmm2,238 ;11101110b-> t,t,z_s,z_r            psrldq xmm2,8        ;SSE2 shift right by BYTES not BITS! t,t,z_s,z_r
       addps  xmm3,xmm2     ;     (t,t,z_s,z_r) + (t,t,y_s,y_r)  -> t,t,S,R
       unpcklps xmm3,xmm3   ;     S,S,R,R  ;!!!!shufps xmm3,xmm3,68  ;     0100 0100b -> S,R,S,R    ;

       X_loop_init:
   mov ecx,[rsp+124]     ;    n_pixX
   movss xmm2,[r10]     ;    *pixelpos[3]
	 movss xmm7,xmm2      ;    t,Z(akt),Y(akt),X(0)

       X_loop:		     ;     !!!2 PIXEL in parallel!!!
	 ;x_p1 & x_p0
	 movss	xmm2,xmm7     ;    pixel  t,t,t,x_p
	 shufps xmm2,xmm2,0   ;    x_p,x_p,x_p,x_p
	 addss	xmm2,xmm4     ;    x_p,x_p,x_p,x_p + t,t,t,resolut -> x_p0,x_p0,x_p0,x_p1
	 shufps xmm2,xmm2,17  ;    0001 0001b ->  x_p1,x_p0,x_p1,x_p0
	 
      ;;;PIPE2
      movaps xmm9,xmm2
      addps  xmm9,xmm10    ;    x_p1,x_p0,x_p1,x_p0+ 2*resolut ->  x_p3,x_p2,x_p3,x_p2 

	 ;x_r &  x_s
	 movss xmm1,xmm6      ;    receiver (trash),z_r,y_r,x_r
	 shufps xmm1,xmm5,0   ;    sender; 0000 0000b -> S,S,R,R

	 subps xmm2,xmm1      ;    x_p-x_r, x_p-x_s -> S1,S0,R1,R0
      ;PIPE2
      subps xmm9,xmm1      ;    x_p-x_r, x_p-x_s -> S3,S2,R3,R2
	 mulps xmm2,xmm2      ;    quadrieren
      ;PIPE2
      mulps xmm9,xmm9      ;    quadrieren
	 
	 addps xmm2,xmm3      ;    S1x,S0x,R1x,R0x + Syz,Syz,Ryz,Ryz
      ;PIPE2
      addps xmm9,xmm3      ;    S1x,S0x,R1x,R0x + Syz,Syz,Ryz,Ryz
	 sqrtps xmm2,xmm2     ;    sqrt
      ;PIPE2
      sqrtps xmm9,xmm9     ;    sqrt

	 movaps xmm1,xmm2     ;
	 shufps xmm1,xmm1,254 ;    1111 1110b -> S1,S1,S1,S0
	 addps xmm2,xmm1     ;    S1,S0,R1,R0 +  S1,S1,S1,S0 = t,t,P1,P0
	 
       ;PIPE2
      movaps xmm11,xmm9     ;
      shufps xmm11,xmm11,254 ;    1111 1110b -> S3,S3,S3,S2
      addps xmm9,xmm11     ;    S3,S2,R3,R2 +  S3,S3,S3,S2 = t,t,P3,P2

	 ;lauflaenge in t zu index
	 mulps xmm2,xmm8      ;    t,t,P1,P0*factor,factor,factor,factor = t,t,index1,index0
	 
      ;PIPE2 lauflaenge in t zu index
      mulps xmm9,xmm8      ;    t,t,P1,P0*factor,factor,factor,factor = t,t,index1,index0

	 ;time to index for P0
	 cvtss2si eax,xmm2    ;    index_value0
   
 	 mov r15,rdi
	 lea rdi,[r13+rbp]
	 movupd xmm14,[r12+rbp] ;read all 2 pixel!
   movupd xmm13,[r12+rbp+16] ;read all 2 pixel! 
   
	;;;;;;;;;;;;;;;;;;;;;;;;
	;rangecheck pixel 1
	 cmp eax,dword [rel WIDTH_BUFFER]	 ;n_Ascan *5      ;cmp eax,dword [rsp+120] ;n_Ascan
	 jge outrange_xsum1	 ;0....2999
	 cmp eax,0		 ;
	 jl outrange_xsum1

	 ;;inrange real
	 ;fld qword [r12+rbp] ; *image_sum_real fpu move
	 shl eax,3		;imul eax,8    ;index-> 64bit (double) addr   NP
	 ;fld qword [r11+rax]	;*buffer_real[index_aktual]         fpu double move
	 ;faddp st1,st0
	 
	 ;;write PIXEL0!
	 ;fstp qword [r13+rbp]	; write pixel to *bild_real[index_aktual]
	 
	 ;in range pixel0
	 ;movss xmm14,[r12+rbp]
	 addsd xmm14,[r11+rax]
	
	   ;complex Ascan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz rangecheck_pixel_2

	   mov rbx,[rsp+40]	 ;*buffer_complex[index_aktual]
	   fld qword [rbx+rax]	 ;fpu move

	   mov rbx,[rsp+104]	 ;*image_sum_compl
	   fld qword [rbx+rbp]
	   faddp st1,st0

	   mov rbx,[rsp+16]	 ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp]  ; rbp = *bild_compl[index_aktual]
	   jmp rangecheck_pixel_2

	outrange_xsum1:
	 ;real
	 fld qword [r12+rbp]  ;*image_sum_real FPU MOVE
	 
	 ;;write pixel0
	 fstp qword [r13+rbp]	 ; write pixel0 to *bild_real[index_aktual]

	   ;Complex AScan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz rangecheck_pixel_2 ;jump out if only real

	   mov rbx,[rsp+104]	 ;*image_sum_compl
	   fld qword [rbx+rbp]

	   mov rbx,[rsp+16]	 ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp]  ; rbp = *bild_compl[index_aktual]

	;;;;;;;;;;;;;;;;;;;;;;;;;;;
	rangecheck_pixel_2:
   
   xor rax,rax          ; clear upper part of RAX for later EAX -> RAX (for delta on memory access)
	 shufps xmm2,xmm2,85	; 0101 0101b -> P1,P1,P1,P1
	 cvtss2si eax,xmm2	; index_value

	 cmp eax,dword [rel WIDTH_BUFFER]	;n_Ascan *5   ;cmp eax,[rsp+120]    ;n_Ascan
	 jge outrange_xsum2	;0....2999
	 cmp eax,0		;
	 jl outrange_xsum2

	 ;;inrange real
	 shl eax,3		;imul eax,8    ;index-> 64bit (double) addr   NP
	 ;fld qword [r12+rbp+8]	;*image_sum_real+ delta Image +1pixel
	 ;fld qword [r11+rax]	;*buffer_real[index_aktual]  fpu double move
	 ;faddp st1,st0
	 
   ;;write pixel1!
	 ;fstp qword [r13+rbp+8] ; write pixel to *bild_real[index_aktual]
   
   ;in range pixel0
	 ;movss xmm14,[r12+rbp]
	 shufpd xmm14,xmm14,1 ; 10   -> rotate
	 addsd xmm14,[r11+rax]
	 shufpd xmm14,xmm14,1 ; 10 -> rotate

	   ;complex Ascan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz out_pixel2

	   mov rbx,[rsp+40]	 ;*buffer_complex[index_aktual]
	   fld qword [rbx+rax]	 ;fpu move

	   mov rbx,[rsp+104]	  ;*image_sum_compl
	   fld qword [rbx+rbp+8]  ; delta Image +1pixel
	   faddp st1,st0

	   mov rbx,[rsp+16]	  ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp+8] ; rbp = *bild_compl[index_aktual]
	   jmp out_pixel2

	outrange_xsum2:
	 ;real
	 fld qword [r12+rbp+8] ;*image_sum_real + deltabild+ 1pixel
	 
	 ;write pixel1!
	 fstp qword [r13+rbp+8] ; write pixel to *bild_real[index_aktual]

	   ;Complex AScan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz out_pixel2	 ;jump out if only real

	   mov rbx,[rsp+104]	  ;*image_sum_compl
	   fld qword [rbx+rbp+8]

	   mov rbx,[rsp+16]	 ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp+8]  ; rbp = *bild_compl[index_aktual]

	;;;;;;;;;;;
	out_pixel2:
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;
	rangecheck_pixel_3:
   
   xor rax,rax          ; clear upper part of RAX for later EAX -> RAX (for delta on memory access)
	 cvtss2si eax,xmm9	; index_value

	 cmp eax,dword [rel WIDTH_BUFFER]	;n_Ascan *5   ;cmp eax,[rsp+120]    ;n_Ascan
	 jge outrange_xsum3	;0....2999
	 cmp eax,0		;
	 jl outrange_xsum3

	 ;inrange real
	 shl eax,3		;imul eax,8    ;index-> 64bit (double) addr   NP
	 ;fld qword [r12+rbp+16]	;*image_sum_real+ delta Image +2pixel
	 ;fld qword [r11+rax]	;*buffer_real[index_aktual]  fpu double move
	 ;faddp st1,st0

	 ;write pixel2!
	 ;fstp qword [r13+rbp+16] ; write pixel to *bild_real[index_aktual]

   ;in range pixel2
	 ;movss xmm14,[r12+rbp]
	 addsd xmm13,[r11+rax]

	   ;complex Ascan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz out_pixel3

	   mov rbx,[rsp+40]	 ;*buffer_complex[index_aktual]
	   fld qword [rbx+rax]	 ;fpu move

	   mov rbx,[rsp+104]	  ;*image_sum_compl
	   fld qword [rbx+rbp+16]  ; delta Image +2pixel
	   faddp st1,st0

	   mov rbx,[rsp+16]	  ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp+16] ; rbp = *bild_compl[index_aktual]
	   jmp out_pixel3

	outrange_xsum3:
	 ;real
	 fld qword [r12+rbp+16] ;*image_sum_real + deltabild+ 2pixel
	 
	 ;write pixel1!
	 fstp qword [r13+rbp+16] ; write pixel to *bild_real[index_aktual]

	   ;Complex AScan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz out_pixel3	 ;jump out if only real

	   mov rbx,[rsp+104]	  ;*image_sum_compl
	   fld qword [rbx+rbp+16]

	   mov rbx,[rsp+16]	 ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp+16]  ; rbp = *bild_compl[index_aktual]

	;;;;;;;;;;;
	out_pixel3:
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;
	rangecheck_pixel_4:
   
   xor rax,rax          ; clear upper part of RAX for later EAX -> RAX (for delta on memory access)
	 shufps xmm9,xmm9,85	; 0101 0101b -> P1,P1,P1,P1
	 cvtss2si eax,xmm9	; index_value

	 cmp eax,dword [rel WIDTH_BUFFER]	;n_Ascan *5   ;cmp eax,[rsp+120]    ;n_Ascan
	 jge outrange_xsum4	;0....2999
	 cmp eax,0		;
	 jl outrange_xsum4

	 ;inrange real
	 shl eax,3		;imul eax,8    ;index-> 64bit (double) addr   NP
	 ;fld qword [r12+rbp+24]	;*image_sum_real+ delta Image +3pixel
	 ;fld qword [r11+rax]	;*buffer_real[index_aktual]  fpu double move
	 ;faddp st1,st0

	 ;;write pixel3!
	 ;fstp qword [r13+rbp+24] ; write pixel to *bild_real[index_aktual]

   ;in range pixel3
	 ;movss xmm14,[r12+rbp]
	 shufpd xmm13,xmm13,1 ; 01 -> rotate
	 addsd xmm13,[r11+rax]
	 shufpd xmm13,xmm13,1 ; 01 -> rotate
	 
	   ;complex Ascan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz out_pixel4

	   mov rbx,[rsp+40]	 ;*buffer_complex[index_aktual]
	   fld qword [rbx+rax]	 ;fpu move

	   mov rbx,[rsp+104]	  ;*image_sum_compl
	   fld qword [rbx+rbp+24]  ; delta Image +3pixel
	   faddp st1,st0

	   mov rbx,[rsp+16]	  ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp+24] ; rbp = *bild_compl[index_aktual]
	   jmp out_pixel4

	outrange_xsum4:
	 ;real
	 fld qword [r12+rbp+24] ;*image_sum_real + deltabild+ 1pixel
	 
	 ;write pixel1!
	 fstp qword [r13+rbp+24] ; write pixel to *bild_real[index_aktual]

	   ;Complex AScan ?
	   cmp rsi,0		 ;dword [rsp+88],0  ;*AScan_complex
	   jz out_pixel4	 ;jump out if only real

	   mov rbx,[rsp+104]	  ;*image_sum_compl
	   fld qword [rbx+rbp+24]

	   mov rbx,[rsp+16]	 ;*bild_compl[index_aktual]
	   fstp qword [rbx+rbp+24]  ; rbp = *bild_compl[index_aktual]

	;;;;;;;;;;;
	out_pixel4:
		
	;real write (NON-temporal)!	
	movupd [rdi],xmm14
	;maskmovdqu xmm14,xmm15
	;lea rdi,[rdi+16]
	movupd [rdi+16],xmm13
	;maskmovdqu xmm13,xmm15
	
	mov rdi,r15	
		

	;add x_p + 4pixel
	addss xmm7,xmm10  ; t,t,t,x_p  + t,t,t,2*resolut
	addss xmm7,xmm10  ; t,t,t,x_p  + t,t,t,2*resolut

	;inc image index (4pixel)
	add rbp,32	 ;*bild_real[index_aktual] +4pixel[double]
	sub  ecx,4
	
check_ifloop:
	cmp  ecx,4
	jge X_loop

	;mod(ecx,4) ==0???
	cmp  ecx,0
	je exit_X_loop  ;finished!!!
	  ;add x size until mod(ecx,4) == 0 -> ecx =4
	  add ecx,1
	  subss xmm7,xmm4      ; t,t,t,x_p  - t,t,t,resolut
	  sub rbp,8	       ;*bild_real[index_aktual] -1pixel[double]
	  jmp check_ifloop

	exit_X_loop:
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;add y_p + 1pixel
      xorps xmm0,xmm0	   ; xmm0 = 0
      movss xmm0,xmm4	   ; t,t,t,resolut
      shufps xmm0,xmm0,81  ; 01010001b -> 0,0,resolut,0
      addps xmm7,xmm0	   ; t,t,y_p,t  + 0,0,resolut,0

      sub  edx,1
      jnz Y_loop
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;add z_p + 1pixel
    xorps xmm0,xmm0	   ; xmm0 = 0
    movss xmm0,xmm4	   ; t,t,t,resolut
    shufps xmm0,xmm0,69    ; 01000101b -> 0,resolut,0,0
    addps xmm7,xmm0	   ; t,z_p,t,t  + 0,resolut,0,0

    sub edi,1
    jnz Z_loop
    jmp exit
    ;;;;end sum_branch;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    exit:
      ;;;;clean up FPU/MMX
	;emms             ;reset mmx, ready fpu
	;FLDCW [esp+4]    ;restore old CTRWORD            ;buffering in stack (16bit!) + FWAIT to ensure save is complete (FSTCW without wait)


      ;;;;;;for debug output
	;mov eax,[rel WIDTH_BUFFER] ;out width (previous test)


	;recover XMM6-7 (callee convention WIN64(?))
	movups xmm7,[rel XMM7_SAVE]
	movups xmm6,[rel XMM6_SAVE]

	;;;;;clean up
	mov rsp, [rel RSP_SAVE]

	;pop eax
	;pop eax

	pop rbx
	pop rdi
	pop rsi
	pop rbp

    
	ret
	;retn ; bug return near!!!


 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;,

	;FNINIT und FLDCW
	;MOVD MM0,eax ; edx nach MM0
	;MOVD MM1,eax ; eax nach MM1
	;PSLLQ MM0,32 ; MM0 um 32 Bit schieben
	;POR MM0,MM1 ; MM1 dazu 'ORen'



      ;inrange:   NOT 2 pipe version
	;sub eax,1          ;1 addr-> 0 addr                         U
	;shl eax,3          ;imul eax,8    ;index-> 64bit addr       U
	;movd ebx,mm1       ;*buffer
	;add eax,ebx        ;*buffer[index]

	;64 bit mmx variante
	;movq mm0, qword [eax]       ;                          U
	;MOVntq qword [edx], mm0     ; writethrough             U
				     ; MMX raised keine flags (I hope)

	;add edx,8      ; + 1 double addr
	;paddd mm2,mm4   ; mm2 + 4  = + 1 int32 addr
	;sub ecx,1
	;jnz fill
	;jmp exit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;32-bit variante
	;mov ebx,dword [eax]     ;sum[ecx]<-buffer[index]
	;mov dword [edx],ebx
	;mov ebx,dword [eax+4]   ;sum[ecx]<-buffer[index]
	;mov dword [edx+4],ebx

      ;64bit mmx memcopy variante
       ;movq mm0, qword [eax]     ;movq nur upipe(?)           U
       ;MOVntq qword [edx], mm0 ;writethrough moveq