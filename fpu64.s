	.file	"fpu.c"
	.text
.globl fpu_53bits
	.type	fpu_53bits, @function
fpu_53bits:
.LFB2:
	movw	$639, -2(%rsp)
#APP
	fldcw -2(%rsp)
#NO_APP
	ret
.LFE2:
	.size	fpu_53bits, .-fpu_53bits
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x1
	.byte	0x3
	.byte	0xc
	.uleb128 0x7
	.uleb128 0x8
	.byte	0x90
	.uleb128 0x1
	.align 8
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB2
	.long	.LFE2-.LFB2
	.uleb128 0x0
	.align 8
.LEFDE1:
	.ident	"GCC: (GNU) 4.2.2 20070909 (prerelease) (4.2.2-0.RC.1mdv2008.0)"
	.section	.note.GNU-stack,"",@progbits
