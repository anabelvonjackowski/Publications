# Install Homebrew
/bin/bash -c "$(curl -fsSL https://githubusercontent.com)"

# Connect Homebrew to Path
(echo; echo 'eval "$(/opt/homebrew/bin/brew shellenv)"') >> ~/.bash_profile
eval "$(/opt/homebrew/bin/brew shellenv)"

# Install FFmepeg
brew install ffmpeg

# verify FFmpeg version
ffmpeg -version

# Extract gif every 0.5 seconds
ffmpeg -i input.mp4 -vf "fps=2" output_%04d.png
