import os, pexpect

def rsync(cmd, pw_file=None):
    #rsync -z --no-R --inplace --files-from=/gpcr/users/mlopez/allosteric/AlloViz_to_GPCRmd/AlloViz/run/slurm_files/dyn9/dyn9.files ori:/protwis/sites/files .

    if not pw_file:
        pw_file = f"{os.path.dirname(os.path.realpath(__file__))}/password.txt" # /gpcr/users/frann/AlloViz_GPCRmd
    
    # print('pw_file',pw_file)
    #/gpcr/users/mlopez/allosteric/AlloViz_to_GPCRmd/AlloViz/run/password.txt
    
    with open(pw_file, "r") as pw:
        child = pexpect.spawn(f'bash -c "{cmd}"', encoding='utf-8')
        # print(child)
        option = child.expect([".* password:",".*Are you sure you want to continue connecting.*"])
        if option == 1:
            child.sendline("yes")
            child.expect(".* password:")
        
        child.sendline(pw.read().strip())
        child.expect(pexpect.EOF, timeout=None)
        out = child.before
        if "error" in out:
            raise Exception(out)
    
    return

