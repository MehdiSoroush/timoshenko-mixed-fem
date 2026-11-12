# Publishing the repository

The local package is a release candidate, not an already published GitHub project.

1. Review the `MIT` license and author attribution in `LICENSE`.
2. Create an empty GitHub repository named `timoshenko-mixed-fem`, without auto-generating a README, gitignore, or license.
3. In PowerShell, open the extracted repository directory and run:

   ```powershell
   git init
   git add .
   git commit -m "Release verified Timoshenko mixed FEM reconstruction v0.4.0"
   git branch -M main
   git remote add origin https://github.com/<YOUR_USERNAME>/timoshenko-mixed-fem.git
   git push -u origin main
   ```

4. Replace `<YOUR_USERNAME>` with your actual GitHub username; do not include the angle brackets.
5. Open the Actions tab and verify that the Python tests pass on both Windows and Ubuntu. Check the Markdown technical report and PDF render correctly.
6. Add the topics listed in [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) and pin this repository to your GitHub profile.
7. Only then add the live GitHub URL to your CV.

**Public claim:** verified Python reconstruction based on a 2017 computational mechanics study, with new 2026 independent checks. **Do not claim:** exact reproduction of the archived Abaqus results.
