	.file	"fpu.c"
	.text
.globl fpu_53bits
	.type	fpu_53bits, @function
fpu_53bits:
	subl	$16, %esp
	movw	$639, 14(%esp)
#APP
	fldcw 14(%esp)
#NO_APP
	addl	$16, %esp
	ret
	.size	fpu_53bits, .-fpu_53bits
	.ident	"GCC: (GNU) 4.2.2 20070909 (prerelease) (4.2.2-0.RC.1mdv2008.0)"
	.section	.note.GNU-stack,"",@progbits
