source "https://rubygems.org"

# --- Core Jekyll Gems

gem "jekyll"

# Jekyll theme
gem "minima", "~> 2.5"

# If you have any plugins, put them here!
group :jekyll_plugins do
  gem "jekyll-scholar"
#  gem "jekyll-autoprefixer", "~> 1.0"
#  gem "jekyll-feed", "~> 0.12"
end

# Windows and JRuby does not include zoneinfo files, so bundle the tzinfo-data gem
# and associated library.
platforms :mingw, :x64_mingw, :mswin, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
end

# Performance-booster for watching directories on Windows
gem "wdm", "~> 0.1.1", :platforms => [:mingw, :x64_mingw, :mswin]

# Lock `http_parser.rb` gem to `v0.6.x` on JRuby builds since newer versions of the gem
# do not have a Java counterpart.
gem "http_parser.rb", "~> 0.6.0", :platforms => [:jruby]

# --- Additional Gems

gem "base64", "~> 0.2.0"

gem "csv", "~> 3.3"

gem "kramdown", "~> 2.4"

gem "logger", "~> 1.6"

gem "observer", "~> 0.1.2"
