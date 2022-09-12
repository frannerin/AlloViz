import os, pexpect

def rsync(cmd, pw_file=None):
    if not pw_file:
        pw_file = f"{os.path.dirname(os.path.realpath(__file__))}/password.txt" # /gpcr/users/frann/AlloViz_GPCRmd
    
    with open(pw_file, "r") as pw:
        child = pexpect.spawn(cmd, encoding='utf-8')
        option = child.expect([".* password:",
                               ".*Are you sure you want to continue connecting.*"
                              ])
        if option == 1:
            child.sendline("yes")
            child.expect(".* password:")
        
        child.sendline(pw.read().strip())
        child.expect(pexpect.EOF, timeout=None)
        out = child.before
        if "error" in out:
            raise Exception(out)
    
    return

